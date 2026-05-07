import { openDB, IDBPDatabase } from 'idb';

export const DEFAULT_CATEGORIES = [
  'Estructura',
  'Instalaciones',
  'Terminaciones',
  'Patologías',
  'Otros',
];

export type Obra = {
  id: string;
  name: string;
  description?: string;
  createdAt: number;
  updatedAt: number;
};

export type Photo = {
  id: string;
  obraId: string;
  blob: Blob;
  thumbBlob?: Blob;
  note: string;
  category: string;
  takenAt: number;
  createdAt: number;
};

export type Settings = {
  id: 'settings';
  categories: string[];
  author?: string;
};

const DB_NAME = 'obras-db';
const DB_VERSION = 1;

let dbPromise: Promise<IDBPDatabase> | null = null;

function getDB() {
  if (!dbPromise) {
    dbPromise = openDB(DB_NAME, DB_VERSION, {
      upgrade(db) {
        if (!db.objectStoreNames.contains('obras')) {
          const obras = db.createObjectStore('obras', { keyPath: 'id' });
          obras.createIndex('updatedAt', 'updatedAt');
        }
        if (!db.objectStoreNames.contains('photos')) {
          const photos = db.createObjectStore('photos', { keyPath: 'id' });
          photos.createIndex('obraId', 'obraId');
          photos.createIndex('createdAt', 'createdAt');
        }
        if (!db.objectStoreNames.contains('settings')) {
          db.createObjectStore('settings', { keyPath: 'id' });
        }
      },
    });
  }
  return dbPromise;
}

function uid(): string {
  if (typeof crypto !== 'undefined' && 'randomUUID' in crypto) {
    return crypto.randomUUID();
  }
  return Math.random().toString(36).slice(2) + Date.now().toString(36);
}

// Settings
export async function getSettings(): Promise<Settings> {
  const db = await getDB();
  const existing = (await db.get('settings', 'settings')) as Settings | undefined;
  if (existing) return existing;
  const initial: Settings = { id: 'settings', categories: [...DEFAULT_CATEGORIES] };
  await db.put('settings', initial);
  return initial;
}

export async function updateSettings(patch: Partial<Settings>) {
  const db = await getDB();
  const current = await getSettings();
  const next: Settings = { ...current, ...patch, id: 'settings' };
  await db.put('settings', next);
  return next;
}

// Obras
export async function listObras(): Promise<Obra[]> {
  const db = await getDB();
  const all = (await db.getAll('obras')) as Obra[];
  return all.sort((a, b) => b.updatedAt - a.updatedAt);
}

export async function getObra(id: string): Promise<Obra | undefined> {
  const db = await getDB();
  return (await db.get('obras', id)) as Obra | undefined;
}

export async function createObra(input: { name: string; description?: string }): Promise<Obra> {
  const db = await getDB();
  const now = Date.now();
  const obra: Obra = {
    id: uid(),
    name: input.name.trim(),
    description: input.description?.trim() || '',
    createdAt: now,
    updatedAt: now,
  };
  await db.put('obras', obra);
  return obra;
}

export async function updateObra(id: string, patch: Partial<Omit<Obra, 'id' | 'createdAt'>>) {
  const db = await getDB();
  const current = (await db.get('obras', id)) as Obra | undefined;
  if (!current) throw new Error('Obra no encontrada');
  const next: Obra = { ...current, ...patch, updatedAt: Date.now() };
  await db.put('obras', next);
  return next;
}

export async function deleteObra(id: string) {
  const db = await getDB();
  const tx = db.transaction(['obras', 'photos'], 'readwrite');
  await tx.objectStore('obras').delete(id);
  const photoStore = tx.objectStore('photos');
  const idx = photoStore.index('obraId');
  let cursor = await idx.openCursor(id);
  while (cursor) {
    await cursor.delete();
    cursor = await cursor.continue();
  }
  await tx.done;
}

// Photos
export async function listPhotosByObra(obraId: string): Promise<Photo[]> {
  const db = await getDB();
  const all = (await db.getAllFromIndex('photos', 'obraId', obraId)) as Photo[];
  return all.sort((a, b) => a.createdAt - b.createdAt);
}

export async function getPhoto(id: string): Promise<Photo | undefined> {
  const db = await getDB();
  return (await db.get('photos', id)) as Photo | undefined;
}

export async function createPhoto(input: {
  obraId: string;
  blob: Blob;
  thumbBlob?: Blob;
  note: string;
  category: string;
  takenAt?: number;
}): Promise<Photo> {
  const db = await getDB();
  const now = Date.now();
  const photo: Photo = {
    id: uid(),
    obraId: input.obraId,
    blob: input.blob,
    thumbBlob: input.thumbBlob,
    note: input.note.trim(),
    category: input.category,
    takenAt: input.takenAt ?? now,
    createdAt: now,
  };
  await db.put('photos', photo);
  await touchObra(input.obraId);
  return photo;
}

export async function updatePhoto(
  id: string,
  patch: Partial<Pick<Photo, 'note' | 'category'>>,
) {
  const db = await getDB();
  const current = (await db.get('photos', id)) as Photo | undefined;
  if (!current) throw new Error('Foto no encontrada');
  const next: Photo = { ...current, ...patch };
  await db.put('photos', next);
  await touchObra(current.obraId);
  return next;
}

export async function deletePhoto(id: string) {
  const db = await getDB();
  const photo = (await db.get('photos', id)) as Photo | undefined;
  if (!photo) return;
  await db.delete('photos', id);
  await touchObra(photo.obraId);
}

async function touchObra(id: string) {
  const db = await getDB();
  const obra = (await db.get('obras', id)) as Obra | undefined;
  if (obra) {
    obra.updatedAt = Date.now();
    await db.put('obras', obra);
  }
}

// Backup
export async function exportAll(): Promise<{ obras: Obra[]; photos: any[]; settings: Settings }> {
  const db = await getDB();
  const obras = (await db.getAll('obras')) as Obra[];
  const photosRaw = (await db.getAll('photos')) as Photo[];
  const settings = await getSettings();
  const photos = await Promise.all(
    photosRaw.map(async (p) => ({
      ...p,
      blob: await blobToDataUrl(p.blob),
      thumbBlob: p.thumbBlob ? await blobToDataUrl(p.thumbBlob) : undefined,
    })),
  );
  return { obras, photos, settings };
}

export async function importAll(payload: {
  obras: Obra[];
  photos: any[];
  settings?: Settings;
}, mode: 'merge' | 'replace' = 'merge') {
  const db = await getDB();
  const tx = db.transaction(['obras', 'photos', 'settings'], 'readwrite');
  if (mode === 'replace') {
    await tx.objectStore('obras').clear();
    await tx.objectStore('photos').clear();
  }
  for (const o of payload.obras || []) {
    await tx.objectStore('obras').put(o);
  }
  for (const p of payload.photos || []) {
    const blob = await dataUrlToBlob(p.blob);
    const thumbBlob = p.thumbBlob ? await dataUrlToBlob(p.thumbBlob) : undefined;
    await tx.objectStore('photos').put({ ...p, blob, thumbBlob });
  }
  if (payload.settings) {
    await tx.objectStore('settings').put({ ...payload.settings, id: 'settings' });
  }
  await tx.done;
}

function blobToDataUrl(blob: Blob): Promise<string> {
  return new Promise((resolve, reject) => {
    const reader = new FileReader();
    reader.onload = () => resolve(reader.result as string);
    reader.onerror = () => reject(reader.error);
    reader.readAsDataURL(blob);
  });
}

async function dataUrlToBlob(dataUrl: string): Promise<Blob> {
  const res = await fetch(dataUrl);
  return await res.blob();
}
