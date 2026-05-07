import jsPDF from 'jspdf';
import type { Obra, Photo } from './db';
import { blobToDataUrl } from './image';

type ReportInput = {
  obra: Obra;
  photos: Photo[];
  reportDescription?: string;
  author?: string;
};

export async function generateReportPdf({
  obra,
  photos,
  reportDescription,
  author,
}: ReportInput): Promise<Blob> {
  const doc = new jsPDF({ unit: 'pt', format: 'a4' });
  const pageW = doc.internal.pageSize.getWidth();
  const pageH = doc.internal.pageSize.getHeight();
  const margin = 40;
  const contentW = pageW - margin * 2;

  drawHeader(doc, obra, author, margin, pageW);
  let y = 110;

  doc.setFont('helvetica', 'bold');
  doc.setFontSize(11);
  doc.setTextColor(60);
  doc.text('Descripción del informe', margin, y);
  y += 14;
  doc.setFont('helvetica', 'normal');
  doc.setFontSize(11);
  doc.setTextColor(20);
  const desc = (reportDescription || '').trim() || '—';
  const descLines = doc.splitTextToSize(desc, contentW);
  doc.text(descLines, margin, y);
  y += descLines.length * 14 + 10;

  if (obra.description) {
    doc.setFont('helvetica', 'bold');
    doc.setFontSize(11);
    doc.setTextColor(60);
    doc.text('Notas de la obra', margin, y);
    y += 14;
    doc.setFont('helvetica', 'normal');
    doc.setTextColor(20);
    const ob = doc.splitTextToSize(obra.description, contentW);
    doc.text(ob, margin, y);
    y += ob.length * 14 + 10;
  }

  doc.setDrawColor(220);
  doc.line(margin, y, pageW - margin, y);
  y += 14;

  // Photos: max two per page in portrait, one large per page if very tall.
  const maxImgW = contentW;
  const maxImgH = 320;

  for (let i = 0; i < photos.length; i++) {
    const p = photos[i];
    const dataUrl = await blobToDataUrl(p.blob);
    const dims = await getImageDims(dataUrl);
    const scale = Math.min(maxImgW / dims.w, maxImgH / dims.h, 1);
    const w = dims.w * scale;
    const h = dims.h * scale;

    const noteLines = doc.splitTextToSize(p.note || '—', contentW);
    const blockH = h + 14 + 14 + noteLines.length * 13 + 18;

    if (y + blockH > pageH - margin) {
      doc.addPage();
      drawPageHeader(doc, obra, margin, pageW);
      y = 70;
    }

    const x = margin + (contentW - w) / 2;
    try {
      doc.addImage(dataUrl, 'JPEG', x, y, w, h, undefined, 'FAST');
    } catch {
      doc.addImage(dataUrl, 'PNG', x, y, w, h, undefined, 'FAST');
    }
    y += h + 8;

    doc.setFont('helvetica', 'bold');
    doc.setFontSize(10);
    doc.setTextColor(60);
    const taken = new Date(p.takenAt).toLocaleString('es-AR');
    doc.text(`Foto ${i + 1} · ${p.category} · ${taken}`, margin, y);
    y += 12;

    doc.setFont('helvetica', 'normal');
    doc.setFontSize(11);
    doc.setTextColor(20);
    doc.text(noteLines, margin, y);
    y += noteLines.length * 13 + 18;
  }

  drawFooters(doc);
  return doc.output('blob');
}

function drawHeader(doc: jsPDF, obra: Obra, author: string | undefined, margin: number, pageW: number) {
  doc.setFillColor(245, 158, 11);
  doc.rect(0, 0, pageW, 8, 'F');

  doc.setFont('helvetica', 'bold');
  doc.setFontSize(18);
  doc.setTextColor(20);
  doc.text('Informe de obra', margin, 38);

  doc.setFont('helvetica', 'bold');
  doc.setFontSize(13);
  doc.text(obra.name, margin, 60);

  doc.setFont('helvetica', 'normal');
  doc.setFontSize(10);
  doc.setTextColor(90);
  const today = new Date().toLocaleDateString('es-AR', {
    year: 'numeric',
    month: 'long',
    day: 'numeric',
  });
  doc.text(`Fecha: ${today}`, margin, 78);
  if (author) doc.text(`Autor: ${author}`, margin, 92);
}

function drawPageHeader(doc: jsPDF, obra: Obra, margin: number, pageW: number) {
  doc.setFillColor(245, 158, 11);
  doc.rect(0, 0, pageW, 4, 'F');
  doc.setFont('helvetica', 'bold');
  doc.setFontSize(11);
  doc.setTextColor(60);
  doc.text(`Informe — ${obra.name}`, margin, 28);
  doc.setDrawColor(230);
  doc.line(margin, 36, pageW - margin, 36);
}

function drawFooters(doc: jsPDF) {
  const pageCount = doc.getNumberOfPages();
  for (let i = 1; i <= pageCount; i++) {
    doc.setPage(i);
    const pageW = doc.internal.pageSize.getWidth();
    const pageH = doc.internal.pageSize.getHeight();
    doc.setFont('helvetica', 'normal');
    doc.setFontSize(9);
    doc.setTextColor(140);
    doc.text(`Página ${i} de ${pageCount}`, pageW / 2, pageH - 18, { align: 'center' });
  }
}

function getImageDims(dataUrl: string): Promise<{ w: number; h: number }> {
  return new Promise((resolve, reject) => {
    const img = new Image();
    img.onload = () => resolve({ w: img.naturalWidth, h: img.naturalHeight });
    img.onerror = (e) => reject(e);
    img.src = dataUrl;
  });
}
