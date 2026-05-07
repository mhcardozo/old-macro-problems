import { useEffect, useRef, useState } from 'react';
import { Link } from 'react-router-dom';
import {
  Obra,
  createObra,
  exportAll,
  importAll,
  listObras,
  getSettings,
  updateSettings,
} from '../lib/db';
import { shareOrDownload } from '../lib/download';

export default function ObrasList() {
  const [obras, setObras] = useState<Obra[]>([]);
  const [showNew, setShowNew] = useState(false);
  const [showSettings, setShowSettings] = useState(false);
  const [toast, setToast] = useState<string | null>(null);
  const fileInputRef = useRef<HTMLInputElement>(null);

  async function refresh() {
    setObras(await listObras());
  }

  useEffect(() => {
    refresh();
  }, []);

  function flash(msg: string) {
    setToast(msg);
    setTimeout(() => setToast(null), 2000);
  }

  async function handleExport() {
    const data = await exportAll();
    const blob = new Blob([JSON.stringify(data)], { type: 'application/json' });
    const stamp = new Date().toISOString().slice(0, 10);
    await shareOrDownload(blob, `obras-backup-${stamp}.json`, 'Backup Obras');
    flash('Backup generado');
  }

  async function handleImport(file: File) {
    try {
      const text = await file.text();
      const json = JSON.parse(text);
      await importAll(json, 'merge');
      await refresh();
      flash('Backup importado');
    } catch (e) {
      console.error(e);
      flash('Error al importar');
    }
  }

  return (
    <div>
      <div className="row" style={{ marginBottom: 12 }}>
        <button className="btn" onClick={() => setShowSettings(true)}>⚙ Ajustes</button>
        <button className="btn" onClick={handleExport}>⬇ Backup</button>
        <button className="btn" onClick={() => fileInputRef.current?.click()}>
          ⬆ Importar
        </button>
        <input
          ref={fileInputRef}
          type="file"
          accept="application/json"
          style={{ display: 'none' }}
          onChange={(e) => {
            const f = e.target.files?.[0];
            if (f) handleImport(f);
            e.target.value = '';
          }}
        />
      </div>

      <h2 className="section-title">Tus obras</h2>

      {obras.length === 0 ? (
        <div className="empty">
          <p>Todavía no creaste ninguna obra.</p>
          <p>Tocá "Nueva obra" para empezar.</p>
        </div>
      ) : (
        <div className="list">
          {obras.map((o) => (
            <Link key={o.id} to={`/obras/${o.id}`} className="card card-link">
              <div>
                <div className="card-title">{o.name}</div>
                <div className="card-sub">
                  Actualizado {new Date(o.updatedAt).toLocaleDateString('es-AR')}
                </div>
              </div>
              <span style={{ color: 'var(--muted)' }}>›</span>
            </Link>
          ))}
        </div>
      )}

      <div className="fab-bar">
        <button className="btn btn-primary btn-block" onClick={() => setShowNew(true)}>
          + Nueva obra
        </button>
      </div>

      {showNew && (
        <NewObraDialog
          onClose={() => setShowNew(false)}
          onCreated={async () => {
            setShowNew(false);
            await refresh();
          }}
        />
      )}

      {showSettings && (
        <SettingsDialog onClose={() => setShowSettings(false)} onSaved={() => flash('Ajustes guardados')} />
      )}

      {toast && <div className="toast">{toast}</div>}
    </div>
  );
}

function NewObraDialog({ onClose, onCreated }: { onClose: () => void; onCreated: () => void }) {
  const [name, setName] = useState('');
  const [description, setDescription] = useState('');
  const [busy, setBusy] = useState(false);

  async function submit() {
    if (!name.trim()) return;
    setBusy(true);
    await createObra({ name, description });
    setBusy(false);
    onCreated();
  }

  return (
    <div className="dialog-backdrop" onClick={onClose}>
      <div className="dialog" onClick={(e) => e.stopPropagation()}>
        <h2>Nueva obra</h2>
        <div className="field">
          <label>Nombre</label>
          <input
            autoFocus
            value={name}
            onChange={(e) => setName(e.target.value)}
            placeholder="Ej: Casa González"
          />
        </div>
        <div className="field">
          <label>Notas (opcional)</label>
          <textarea
            value={description}
            onChange={(e) => setDescription(e.target.value)}
            placeholder="Dirección, cliente, contacto…"
          />
        </div>
        <div className="row">
          <button className="btn" onClick={onClose}>Cancelar</button>
          <button className="btn btn-primary" disabled={!name.trim() || busy} onClick={submit}>
            Crear
          </button>
        </div>
      </div>
    </div>
  );
}

function SettingsDialog({ onClose, onSaved }: { onClose: () => void; onSaved: () => void }) {
  const [author, setAuthor] = useState('');
  const [categoriesText, setCategoriesText] = useState('');
  const [loaded, setLoaded] = useState(false);

  useEffect(() => {
    getSettings().then((s) => {
      setAuthor(s.author || '');
      setCategoriesText(s.categories.join('\n'));
      setLoaded(true);
    });
  }, []);

  async function save() {
    const categories = categoriesText
      .split('\n')
      .map((s) => s.trim())
      .filter(Boolean);
    await updateSettings({ author: author.trim(), categories });
    onSaved();
    onClose();
  }

  return (
    <div className="dialog-backdrop" onClick={onClose}>
      <div className="dialog" onClick={(e) => e.stopPropagation()}>
        <h2>Ajustes</h2>
        {!loaded ? (
          <p>Cargando…</p>
        ) : (
          <>
            <div className="field">
              <label>Autor (aparece en el PDF)</label>
              <input value={author} onChange={(e) => setAuthor(e.target.value)} placeholder="Tu nombre" />
            </div>
            <div className="field">
              <label>Categorías (una por línea)</label>
              <textarea
                value={categoriesText}
                onChange={(e) => setCategoriesText(e.target.value)}
                style={{ minHeight: 140 }}
              />
            </div>
            <div className="row">
              <button className="btn" onClick={onClose}>Cancelar</button>
              <button className="btn btn-primary" onClick={save}>Guardar</button>
            </div>
          </>
        )}
      </div>
    </div>
  );
}
