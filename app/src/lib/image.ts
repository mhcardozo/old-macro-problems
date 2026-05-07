// Resize an image File/Blob to fit within maxDim px keeping aspect ratio.
// Returns a JPEG Blob plus its data URL for previewing.
export async function processImage(
  file: File | Blob,
  maxDim = 1600,
  quality = 0.82,
): Promise<{ blob: Blob; thumb: Blob; width: number; height: number }> {
  const bitmap = await createImageBitmapFromBlob(file);
  const { width, height } = fit(bitmap.width, bitmap.height, maxDim);
  const blob = await drawToJpeg(bitmap, width, height, quality);
  const thumbSize = fit(bitmap.width, bitmap.height, 480);
  const thumb = await drawToJpeg(bitmap, thumbSize.width, thumbSize.height, 0.7);
  bitmap.close?.();
  return { blob, thumb, width, height };
}

export function blobToObjectUrl(blob: Blob): string {
  return URL.createObjectURL(blob);
}

async function createImageBitmapFromBlob(file: File | Blob): Promise<ImageBitmap> {
  // imageOrientation: 'from-image' applies EXIF rotation when supported.
  try {
    return await createImageBitmap(file, { imageOrientation: 'from-image' } as any);
  } catch {
    return await createImageBitmap(file);
  }
}

function fit(w: number, h: number, max: number) {
  if (w <= max && h <= max) return { width: w, height: h };
  const r = w > h ? max / w : max / h;
  return { width: Math.round(w * r), height: Math.round(h * r) };
}

async function drawToJpeg(
  bitmap: ImageBitmap,
  w: number,
  h: number,
  quality: number,
): Promise<Blob> {
  const canvas = document.createElement('canvas');
  canvas.width = w;
  canvas.height = h;
  const ctx = canvas.getContext('2d')!;
  ctx.drawImage(bitmap, 0, 0, w, h);
  return await new Promise<Blob>((resolve) =>
    canvas.toBlob((b) => resolve(b!), 'image/jpeg', quality),
  );
}

export function blobToDataUrl(blob: Blob): Promise<string> {
  return new Promise((resolve, reject) => {
    const reader = new FileReader();
    reader.onload = () => resolve(reader.result as string);
    reader.onerror = () => reject(reader.error);
    reader.readAsDataURL(blob);
  });
}
