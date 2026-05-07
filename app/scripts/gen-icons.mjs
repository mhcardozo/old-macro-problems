// Generate simple PNG icons (192 and 512) without external deps.
// Pure-JS minimal PNG encoder (8-bit RGBA, no compression — uses zlib stored blocks).
import { writeFileSync } from 'node:fs';
import { dirname, join } from 'node:path';
import { fileURLToPath } from 'node:url';
import zlib from 'node:zlib';

const __dirname = dirname(fileURLToPath(import.meta.url));
const outDir = join(__dirname, '..', 'public');

function buildIcon(size) {
  const data = Buffer.alloc(size * size * 4);
  const bg = [245, 158, 11, 255]; // accent orange
  const fg = [26, 17, 0, 255];     // dark text

  // Fill background.
  for (let i = 0; i < size * size; i++) {
    data[i * 4] = bg[0];
    data[i * 4 + 1] = bg[1];
    data[i * 4 + 2] = bg[2];
    data[i * 4 + 3] = bg[3];
  }

  // Draw a thick "O" (ring) centered.
  const cx = size / 2;
  const cy = size / 2;
  const rOuter = size * 0.34;
  const rInner = size * 0.22;
  for (let y = 0; y < size; y++) {
    for (let x = 0; x < size; x++) {
      const dx = x - cx;
      const dy = y - cy;
      const d = Math.sqrt(dx * dx + dy * dy);
      if (d <= rOuter && d >= rInner) {
        const i = (y * size + x) * 4;
        data[i] = fg[0];
        data[i + 1] = fg[1];
        data[i + 2] = fg[2];
        data[i + 3] = fg[3];
      }
    }
  }

  return encodePNG(size, size, data);
}

function encodePNG(width, height, rgba) {
  // PNG signature
  const sig = Buffer.from([137, 80, 78, 71, 13, 10, 26, 10]);

  // IHDR
  const ihdr = Buffer.alloc(13);
  ihdr.writeUInt32BE(width, 0);
  ihdr.writeUInt32BE(height, 4);
  ihdr[8] = 8;  // bit depth
  ihdr[9] = 6;  // color type RGBA
  ihdr[10] = 0; // compression
  ihdr[11] = 0; // filter
  ihdr[12] = 0; // interlace

  // IDAT: scanlines with filter byte 0
  const stride = width * 4;
  const raw = Buffer.alloc((stride + 1) * height);
  for (let y = 0; y < height; y++) {
    raw[y * (stride + 1)] = 0;
    rgba.copy(raw, y * (stride + 1) + 1, y * stride, y * stride + stride);
  }
  const compressed = zlib.deflateSync(raw);

  // IEND
  return Buffer.concat([sig, chunk('IHDR', ihdr), chunk('IDAT', compressed), chunk('IEND', Buffer.alloc(0))]);
}

function chunk(type, data) {
  const len = Buffer.alloc(4);
  len.writeUInt32BE(data.length, 0);
  const typeBuf = Buffer.from(type, 'ascii');
  const crc = Buffer.alloc(4);
  crc.writeUInt32BE(crc32(Buffer.concat([typeBuf, data])) >>> 0, 0);
  return Buffer.concat([len, typeBuf, data, crc]);
}

const crcTable = (() => {
  const t = new Uint32Array(256);
  for (let n = 0; n < 256; n++) {
    let c = n;
    for (let k = 0; k < 8; k++) c = c & 1 ? 0xedb88320 ^ (c >>> 1) : c >>> 1;
    t[n] = c >>> 0;
  }
  return t;
})();
function crc32(buf) {
  let c = 0xffffffff;
  for (let i = 0; i < buf.length; i++) c = crcTable[(c ^ buf[i]) & 0xff] ^ (c >>> 8);
  return (c ^ 0xffffffff) >>> 0;
}

writeFileSync(join(outDir, 'icon-192.png'), buildIcon(192));
writeFileSync(join(outDir, 'icon-512.png'), buildIcon(512));
console.log('Generated icon-192.png and icon-512.png');
