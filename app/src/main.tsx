import React from 'react';
import ReactDOM from 'react-dom/client';
import { HashRouter, Route, Routes, Navigate } from 'react-router-dom';
import App from './App';
import ObrasList from './pages/ObrasList';
import ObraDetail from './pages/ObraDetail';
import PhotoEdit from './pages/PhotoEdit';
import './index.css';

ReactDOM.createRoot(document.getElementById('root')!).render(
  <React.StrictMode>
    <HashRouter>
      <Routes>
        <Route path="/" element={<App />}>
          <Route index element={<ObrasList />} />
          <Route path="obras/:obraId" element={<ObraDetail />} />
          <Route path="obras/:obraId/photos/new" element={<PhotoEdit />} />
          <Route path="obras/:obraId/photos/:photoId" element={<PhotoEdit />} />
          <Route path="*" element={<Navigate to="/" replace />} />
        </Route>
      </Routes>
    </HashRouter>
  </React.StrictMode>,
);
