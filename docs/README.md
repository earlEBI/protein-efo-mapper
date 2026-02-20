# GitHub Pages for pQTL Mapper

This folder is a static GitHub Pages site.

## Publish

1. Push this repository to GitHub.
2. In GitHub, open `Settings -> Pages`.
3. Ensure GitHub Pages is enabled.
4. The workflow `.github/workflows/pages.yml` will deploy this `docs/` folder.

After deployment, your page URL will be:

`https://<your-user-or-org>.github.io/<repo-name>/`

## Note

GitHub Pages is static hosting only.  
The actual mapping runs in:

- `skills/pqtl-measurement-mapper/scripts/map_measurement_efo.py`
- or your local/API-backed web app.
