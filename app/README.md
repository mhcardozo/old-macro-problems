# Obras — Informes de obra

PWA para tomar fotos en obra, agregar notas y categoría, y generar informes en PDF.

## Características

- Crear varias **obras** (proyectos) y agregar fotos a cada una.
- Tomar foto desde la cámara o elegir de la galería.
- Cada foto lleva **categoría** (editable: Estructura, Instalaciones, Terminaciones, Patologías, Otros…) y **nota**.
- Generar **PDF** del informe con descripción y todas las fotos con sus notas.
- **Backup manual**: exportar todo (obras + fotos) a un JSON, y volver a importarlo.
- Datos guardados localmente en el dispositivo (IndexedDB). Funciona offline.
- Instalable en iPhone (Safari → Compartir → Agregar a pantalla de inicio) y Android (Chrome → Instalar app).

## Desarrollo

```bash
npm install
npm run dev
```

Abre <http://localhost:5173>.

### Probar desde el celular (cámara)

La cámara web requiere HTTPS. Para probar desde un teléfono mientras desarrollás, exponé el dev server con un tunnel:

```bash
# en otra terminal
npm run tunnel
```

Eso usa LocalTunnel y te da una URL HTTPS pública (`https://<algo>.loca.lt`). Abrila desde tu teléfono.

## Build

```bash
npm run build
npm run preview
```

El build genera el service worker y los archivos estáticos en `dist/`.

Para que la PWA quede en una sub-ruta (como GitHub Pages), pasá `BASE_PATH`:

```bash
BASE_PATH=/obras-app/ npm run build
```

## Deploy a GitHub Pages

Hay un workflow en `.github/workflows/deploy-pages.yml` que ya construye y publica
la app cada vez que se pushea a `main` o a la rama de desarrollo.

Para activarlo (una sola vez):

1. En GitHub: **Settings → Pages → Build and deployment → Source: GitHub Actions**.
2. Pushear cualquier cambio dentro de `app/`. El workflow corre y publica.
3. La app queda en `https://<usuario>.github.io/obras-app/`.

Una vez online, abrila desde el iPhone con Safari → **Compartir → Agregar a pantalla de inicio**.

## Datos

- Todo se guarda en IndexedDB del navegador. **No hay servidor.**
- Para mover datos a otro dispositivo, usá Backup ⬇ (genera un `.json`) y luego Importar ⬆ en el nuevo.
- Si borrás los datos del sitio en el navegador, perdés las obras. Hacé backups periódicos.

## Iconos

Se generan con `node scripts/gen-icons.mjs`. Reemplazá `public/icon-192.png` y `public/icon-512.png` por tu propio diseño cuando quieras.
