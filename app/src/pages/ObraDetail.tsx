import { useEffect, useMemo, useRef, useState } from 'react';
import { Link, useNavigate, useParams } from 'react-router-dom';
import {
  Obra,
  Photo,
  deleteObra,
  getObra,
  getSettings,
  listPhotosByObra,
  updateObra,
} from '../lib/db';
import { generateReportPdf } from '../lib/pdf';
import { shareOrDownload } from '../lib/download';
import { blobToObjectUrl } from '../lib/image';

export default function ObraDetail() {
  const { obraId } = useParams<{ obraId: string }>();
  const navigate = useNavigate();
  const [obra, setObra] = useState<Obra | null>(null);
  const [photos, setPhotos] = useState<Photo[]>([]);
  const [showPdfDialog, setShowPdfDialog] = useState(false);
  const [busy, setBusy] = useState(false);
  const [toast, setToast] = useState<string | null>(null);
  const urlsRef = useRef<string[]>([]);

  async function refresh() {
    if (!obraId) return;
    const [o, ps] = await Promise.all([getObra(obraId), listPhotosByObra(obraId)]);
    if (!o) {
      navigate('/', { replace: true });
      return;
    }
    setObra(o);
    setPhotos(ps);
  }

  useEffect(() => {
    refresh();
    return () => {
      urlsRef.current.forEach(URL.revokeObjectURL);
      urlsRef.current = [];
    };
  }, [obraId]);

  const tileUrls = useMemo(() => {
    urlsRef.current.forEach(URL.revokeObjectURL);
    const urls = photos.map((p) => blobToObjectUrl(p.thumbBlob || p.blob));
    urlsRef.current = urls;
    return urls;
  }, [photos]);

  function flash(msg: string) {
    setToast(msg);
    setTimeout(() => setToast(null), 2000);
  }

  async function handleDelete() {
    if (!obra) return;
    if (!confirm(`¿Eliminar la obra "${obra.name}" y todas sus fotos?`)) return;
    await deleteObra(obra.id);
    navigate('/', { replace: true });
  }

  async function handleRename() {
    if (!obra) return;
    const next = prompt('Nuevo nombre de la obra:', obra.name);
    if (next === null) return;
    const trimmed = next.trim();
    if (!trimmed) return;
    await updateObra(obra.id, { name: trimmed });
    refresh();
  }

  async function handleEditNotes() {
    if (!obra) return;
    const next = prompt('Notas de la obra:', obra.description || '');
    if (next === null) return;
    await updateObra(obra.id, { description: next });
    refresh();
  }

  if (!obra) return <div className="empty">Cargando…</div>;

  return (
    <div>
      <div className="card">
        <div className="card-title" style={{ fontSize: 18 }}>{obra.name}</div>
        {obra.description && <div className="card-sub" style={{ marginTop: 6 }}>{obra.description}</div>}
        <div className="row" style={{ marginTop: 12 }}>
          <button className="btn" onClick={handleRename}>Renombrar</button>
          <button className="btn" onClick={handleEditNotes}>Notas</button>
          <button className="btn btn-danger" onClick={handleDelete}>Eliminar</button>
        </div>
      </div>

      <h2 className="section-title">Fotos ({photos.length})</h2>

      {photos.length === 0 ? (
        <div className="empty">
          <p>No hay fotos todavía.</p>
        </div>
      ) : (
        <div className="photo-grid">
          {photos.map((p, i) => (
            <Link key={p.id} to={`/obras/${obra.id}/photos/${p.id}`} className="photo-tile">
              <img src={tileUrls[i]} alt={p.note || 'Foto'} loading="lazy" />
              <div className="badge">{p.category}</div>
            </Link>
          ))}
        </div>
      )}

      <div className="fab-bar">
        <button
          className="btn btn-primary btn-block"
          onClick={() => navigate(`/obras/${obra.id}/photos/new`)}
        >
          + Tomar / agregar foto
        </button>
        <button
          className="btn btn-block"
          disabled={photos.length === 0}
          onClick={() => setShowPdfDialog(true)}
        >
          Generar informe PDF
        </button>
      </div>

      {showPdfDialog && (
        <PdfDialog
          obra={obra}
          photos={photos}
          busy={busy}
          setBusy={setBusy}
          onClose={() => setShowPdfDialog(false)}
          onDone={(msg) => {
            setShowPdfDialog(false);
            flash(msg);
          }}
        />
      )}

      {toast && <div className="toast">{toast}</div>}
    </div>
  );
}

function PdfDialog({
  obra,
  photos,
  busy,
  setBusy,
  onClose,
  onDone,
}: {
  obra: Obra;
  photos: Photo[];
  busy: boolean;
  setBusy: (b: boolean) => void;
  onClose: () => void;
  onDone: (msg: string) => void;
}) {
  const [reportDescription, setReportDescription] = useState('');
  const [author, setAuthor] = useState('');

  useEffect(() => {
    getSettings().then((s) => setAuthor(s.author || ''));
  }, []);

  async function generate() {
    setBusy(true);
    try {
      const blob = await generateReportPdf({ obra, photos, reportDescription, author });
      const stamp = new Date().toISOString().slice(0, 10);
      const safeName = obra.name.replace(/[^a-zA-Z0-9-_]+/g, '_');
      const result = await shareOrDownload(blob, `informe_${safeName}_${stamp}.pdf`, `Informe ${obra.name}`);
      onDone(result === 'shared' ? 'PDF compartido' : 'PDF generado');
    } catch (e) {
      console.error(e);
      onDone('Error al generar PDF');
    } finally {
      setBusy(false);
    }
  }

  return (
    <div className="dialog-backdrop" onClick={onClose}>
      <div className="dialog" onClick={(e) => e.stopPropagation()}>
        <h2>Generar informe PDF</h2>
        <div className="field">
          <label>Autor</label>
          <input value={author} onChange={(e) => setAuthor(e.target.value)} placeholder="Tu nombre" />
        </div>
        <div className="field">
          <label>Descripción del informe</label>
          <textarea
            value={reportDescription}
            onChange={(e) => setReportDescription(e.target.value)}
            placeholder="Resumen de la visita, observaciones generales…"
          />
        </div>
        <div className="row">
          <button className="btn" onClick={onClose} disabled={busy}>Cancelar</button>
          <button className="btn btn-primary" onClick={generate} disabled={busy}>
            {busy ? 'Generando…' : 'Generar'}
          </button>
        </div>
      </div>
    </div>
  );
}
