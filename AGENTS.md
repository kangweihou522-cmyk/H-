# Repository working instructions

Use the existing checkout; do not create worktrees unless requested. Preserve user inputs, analysis results, local configuration and Docker credentials. Tests and examples must use synthetic data; customer reports and real sequencing data are not published.

## Keep the user guide synchronized

The user explicitly requires complete operating instructions, methods and references in the final workstation navigation item and a matching Word document in GitHub, updated with subsequent changes.

- `docs/USER_GUIDE.md` is the sole editable full guide. When changing analysis methods, thresholds, inputs, outputs, plotting controls, navigation or deployment, review and update the relevant passages and references in this file.
- Generate the matching Word, offline HTML and version manifest together: `python -m pip install -r docs/requirements.txt`, then `python scripts/build_user_guide.py`.
- Commit `docs/USER_GUIDE.md`, `docs/USER_GUIDE.docx`, `static/guide/index.html` and `docs/guide-manifest.json` with the code change. Do not manually edit the generated files.
- Before committing, run `python scripts/build_user_guide.py --check`. GitHub Actions runs the same check and rejects missing/stale artifacts; application source fingerprints deliberately require regeneration when tracked runtime files change. This verifies artifact consistency, while the author must still verify scientific wording against the implementation.
- Preserve the two DEG profiles, clearly distinguish exploratory nominal P from BH FDR, and document skipped or unsupported analyses accurately. Do not describe reference projects as directly executed dependencies unless the code actually uses them.
- The guide must teach beginners: introduce biological/sequencing concepts and define terms before relying on acronyms; retain worked calculations, diagram captions, interpretation examples and exercises. Toy P/padj values must be statistically coherent and explicitly marked as teaching data.
- Teaching image sources are in `docs/assets`; `scripts/build_guide_figures.py` redraws the original synthetic diagrams. The guide builder embeds referenced images in Word and copies them to `static/guide/assets`. Commit source images and generated copies together; do not hand-edit generated HTML/Word.
- Changes on GitHub do not automatically update the user's installed Windows workstation. Tell the user to finish analyses, update via `更新工作站.bat`, and restart the backend to see new code and documentation.

## Focused validation

Use tests appropriate to the change. Guide checks need document-generation dependencies, not the full Docker/R pipeline. Frontend changes should be checked in a real browser, including the manual without any selected project and at mobile width. Do not rerun successful sequencing analyses for a documentation-only change.
