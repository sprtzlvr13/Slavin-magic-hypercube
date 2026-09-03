import struct
import math
import numpy as np

def parse_fits_header(f):
    cards = []
    while True:
        block = f.read(2880)
        for i in range(0, 2880, 80):
            card = block[i:i+80].decode('ascii', errors='replace')
            cards.append(card)
            if card.startswith('END'):
                return cards

def header_dict(cards):
    d = {}
    for c in cards:
        if '=' in c and not c.startswith('COMMENT') and not c.startswith('HISTORY'):
            key = c[:8].strip()
            val = c[9:].split('/')[0].strip()
            d[key] = val
    return d

TFORM_MAP = {'D': ('d', 8), 'E': ('f', 4), 'J': ('i', 4), 'I': ('h', 2), 'K': ('q', 8), 'B': ('B', 1)}

def read_lightcurve_fits(path):
    with open(path, 'rb') as f:
        h0 = parse_fits_header(f)
        d0 = header_dict(h0)
        naxis0 = int(d0.get('NAXIS', 0))
        if naxis0 > 0:
            nbytes = abs(int(d0.get('BITPIX', 8))) // 8
            n = 1
            for i in range(1, naxis0+1):
                n *= int(d0[f'NAXIS{i}'])
            databytes = n * nbytes
            pad = (2880 - databytes % 2880) % 2880
            f.seek(databytes + pad, 1)

        h1 = parse_fits_header(f)
        d1 = header_dict(h1)
        nrows = int(d1['NAXIS2'])
        rowbytes = int(d1['NAXIS1'])
        tfields = int(d1['TFIELDS'])
        cols = []
        offset = 0
        for i in range(1, tfields+1):
            ttype = d1.get(f'TTYPE{i}', '').strip("'").strip()
            tform = d1.get(f'TFORM{i}', '').strip("'").strip()
            j = 0
            while j < len(tform) and tform[j].isdigit():
                j += 1
            repeat = int(tform[:j]) if j > 0 else 1
            typecode = tform[j:j+1]
            cols.append((ttype, typecode, repeat, offset))
            fmt, size = TFORM_MAP.get(typecode, ('B', 1))
            offset += size * repeat

        raw = f.read(nrows * rowbytes)

        def read_col(name):
            for ttype, typecode, repeat, off in cols:
                if ttype == name:
                    fmt, size = TFORM_MAP[typecode]
                    out = np.empty(nrows)
                    for r in range(nrows):
                        base = r*rowbytes + off
                        val = struct.unpack('>' + fmt*repeat, raw[base:base+size*repeat])
                        out[r] = val[0]
                    return out
            raise KeyError(name)

        time = read_col('TIME')
        quality = read_col('QUALITY').astype(int)
        pdcsap = read_col('PDCSAP_FLUX')
    return time, quality, pdcsap

IMPULSIVE_OUTLIER_BIT = 512

if __name__ == "__main__":
    import glob
    paths = {
        'pi_men_s01': '/home/claude/pimen_extract/pi_men_tess_claude_package/original_fits/tess_sector01_pi_men_lc.fits',
        'pi_men_s04': '/home/claude/pimen_extract/pi_men_tess_claude_package/original_fits/tess_sector04_pi_men_lc.fits',
        'pi_men_s08': '/home/claude/pimen_extract/pi_men_tess_claude_package/original_fits/tess_sector08_pi_men_lc.fits',
        'wasp18_s02': '/home/claude/wasp18_extract/wasp18_tess_claude_package/original_fits/tess_sector02_wasp18_lc.fits',
        'wasp18_s03': '/home/claude/wasp18_extract/wasp18_tess_claude_package/original_fits/tess_sector03_wasp18_lc.fits',
    }
    for name, p in paths.items():
        t, q, flux = read_lightcurve_fits(p)
        n_impulsive = int(((q.astype(int) & IMPULSIVE_OUTLIER_BIT) != 0).sum())
        n_any = int((q != 0).sum())
        print(f"{name}: {len(t)} rows, {n_any} any-flag, {n_impulsive} ImpulsiveOutlier(512)")
