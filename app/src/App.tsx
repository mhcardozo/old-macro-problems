import { Outlet, useLocation, useNavigate } from 'react-router-dom';

export default function App() {
  const location = useLocation();
  const navigate = useNavigate();
  const isRoot = location.pathname === '/' || location.pathname === '';

  return (
    <div className="app-shell">
      <header className="app-header">
        {!isRoot && (
          <button className="icon-btn" onClick={() => navigate(-1)} aria-label="Volver">
            ←
          </button>
        )}
        <h1>Obras</h1>
      </header>
      <main className="app-main">
        <Outlet />
      </main>
      <footer className="app-footer">
        <span>v0.1 · datos guardados en este dispositivo</span>
      </footer>
    </div>
  );
}
