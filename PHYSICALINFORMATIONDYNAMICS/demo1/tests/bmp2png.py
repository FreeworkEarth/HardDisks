#!/usr/bin/env python3
"""PID Demo 1c: convert the window's screenshots (24-bit uncompressed BMP, SDL_SaveBMP) to PNG with the standard library only
(struct, zlib), so that no image package is needed on either machine.  Usage: bmp2png.py in.bmp out.png"""
import struct
import sys
import zlib


def read_bmp(path):
    with open(path, "rb") as f:
        data = f.read()
    if data[:2] != b"BM":
        raise SystemExit(f"{path}: not a BMP")
    off, = struct.unpack_from("<I", data, 10)
    hsize, w, h, planes, bpp, comp = struct.unpack_from("<IiiHHI", data, 14)
    if bpp != 24 or comp != 0:
        raise SystemExit(f"{path}: {bpp}-bit, compression {comp}; only 24-bit uncompressed BMP is read")
    top_down = h < 0
    h = abs(h)
    stride = (w * 3 + 3) & ~3
    rows = []
    for y in range(h):
        src = y if top_down else h - 1 - y
        row = data[off + src * stride: off + src * stride + w * 3]
        rgb = bytearray(row)
        rgb[0::3], rgb[2::3] = row[2::3], row[0::3]     # BMP stores B, G, R
        rows.append(bytes(rgb))
    return w, h, rows


def write_png(path, w, h, rows):
    def chunk(kind, body):
        return struct.pack(">I", len(body)) + kind + body + struct.pack(">I", zlib.crc32(kind + body) & 0xFFFFFFFF)
    raw = b"".join(b"\x00" + r for r in rows)          # filter type 0 on every row
    png = b"\x89PNG\r\n\x1a\n" + chunk(b"IHDR", struct.pack(">IIBBBBB", w, h, 8, 2, 0, 0, 0))
    png += chunk(b"IDAT", zlib.compress(raw, 9)) + chunk(b"IEND", b"")
    with open(path, "wb") as f:
        f.write(png)


if __name__ == "__main__":
    if len(sys.argv) != 3:
        raise SystemExit("usage: bmp2png.py in.bmp out.png")
    w, h, rows = read_bmp(sys.argv[1])
    write_png(sys.argv[2], w, h, rows)
    print(f"{sys.argv[2]}: {w}x{h}")
