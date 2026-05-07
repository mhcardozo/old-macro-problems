import { useEffect, useRef, useState } from 'react';
import { useNavigate, useParams } from 'react-router-dom';
import {
  Photo,
  createPhoto,
  deletePhoto,
  getPhoto,
  getSettings,
  updatePhoto,
} from '../lib/db';
import { blobToObjectUrl, processImage } from '../lib/image';

export default function PhotoEdit() {
  const { obraId, photoId } = useParams<{ obraId: string; photoId?: string }>();
  const navigate = useNavigate();
  const isNew = !photoId;

  const [categories, setCategories] = useState<string[]>([]);
  const [category, setCategory] = useState('');
  const [note, setNote] = useState('');
  const [previewUrl, setPreviewUrl] = useState<string | null>(null);
  const [pendingBlobs, setPendingBlobs] = useState<{ blob: Blob; thumb: Blob } | null>(null);
  const [existing, setExisting] = useState<Photo | null>(null);
  const [busy, setBusy] = useState(false);
  const cameraRef = useRef<HTMLInputElement>(null);
  const galleryRef = useRef<HTMLInputElement>(null);
  const objectUrlsRef = useRef<string[]>([]);

  useEffect(() => {
    getSettings().then((s) => {
      setCategories(s.categories);
      if (!category && s.categories.length > 0) setCategory(s.categories[0]);
    });
  }, []);

  useEffect(() => {
    if (!photoId) return;
    getPhoto(photoId).then((p) => {
      if (!p) {
        navigate(`/obras/${obraId}`, { replace: true });
        return;
      }
      setExisting(p);
      setNote(p.note);
      setCategory(p.category);
      const url = blobToObjectUrl(p.blob);
      objectUrlsRef.current.push(url);
      setPreviewUrl(url);
    });
  }, [photoId]);

  useEffect(() => {
    return () => {
      objectUrlsRef.current.forEach(URL.revokeObjectURL);
      objectUrlsRef.current = [];
    };
  }, []);

  async function onFileSelected(file: File | undefined) {
    if (!file) return;
    setBusy(true);
    try {
      const { blob, thumb } = await processImage(file, 1600, 0.82);
      const url = blobToObjectUrl(blob);
      objectUrlsRef.current.push(url);
      setPreviewUrl(url);
      setPendingBlobs({ blob, thumb });
    } catch (e) {
      console.error(e);
      alert('No se pudo procesar la imagen');
    } finally {
      setBusy(false);
    }
  }

  async function save() {
    if (!obraId) return;
    setBusy(true);
    try {
      if (isNew) {
        if (!pendingBlobs) {
          alert('Tomá o elegí una foto primero.');
          setBusy(false);
          return;
        }
        await createPhoto({
          obraId,
          blob: pendingBlobs.blob,
          thumbBlob: pendingBlobs.thumb,
          note,
          category,
        });
      } else if (existing) {
        await updatePhoto(existing.id, { note, category });
      }
      navigate(`/obras/${obraId}`, { replace: true });
    } finally {
      setBusy(false);
    }
  }

  async function handleDelete() {
    if (!existing) return;
    if (!confirm('¿Eliminar esta foto?')) return;
    await deletePhoto(existing.id);
    navigate(`/obras/${obraId}`, { replace: true });
  }

  return (
    <div>
      <div className="preview" style={{ minHeight: 200, display: previewUrl ? 'block' : 'flex', alignItems: 'center', justifyContent: 'center' }}>
        {previewUrl ? (
          <img src={previewUrl} alt="Vista previa" />
        ) : (
          <div style={{ padding: 32, color: 'var(--muted)' }}>Sin imagen</div>
        )}
      </div>

      {isNew && (
        <div className="row" style={{ marginBottom: 12 }}>
          <button className="btn" onClick={() => cameraRef.current?.click()} disabled={busy}>
            📷 Cámara
          </button>
          <button className="btn" onClick={() => galleryRef.current?.click()} disabled={busy}>
            🖼 Galería
          </button>
          <input
            ref={cameraRef}
            type="file"
            accept="image/*"
            capture="environment"
            style={{ display: 'none' }}
            onChange={(e) => {
              onFileSelected(e.target.files?.[0]);
              e.target.value = '';
            }}
          />
          <input
            ref={galleryRef}
            type="file"
            accept="image/*"
            style={{ display: 'none' }}
            onChange={(e) => {
              onFileSelected(e.target.files?.[0]);
              e.target.value = '';
            }}
          />
        </div>
      )}

      <div className="field">
        <label>Categoría</label>
        <select value={category} onChange={(e) => setCategory(e.target.value)}>
          {categories.map((c) => (
            <option key={c} value={c}>
              {c}
            </option>
          ))}
          {category && !categories.includes(category) && (
            <option value={category}>{category}</option>
          )}
        </select>
      </div>

      <div className="field">
        <label>Nota / descripción</label>
        <textarea
          value={note}
          onChange={(e) => setNote(e.target.value)}
          placeholder="Qué muestra la foto, observaciones, mediciones…"
        />
      </div>

      <div className="fab-bar">
        <button
          className="btn btn-primary btn-block"
          onClick={save}
          disabled={busy || (isNew && !pendingBlobs)}
        >
          {busy ? 'Guardando…' : 'Guardar'}
        </button>
        {!isNew && (
          <button className="btn btn-danger btn-block" onClick={handleDelete} disabled={busy}>
            Eliminar foto
          </button>
        )}
      </div>
    </div>
  );
}
