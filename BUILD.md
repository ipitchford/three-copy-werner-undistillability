# Build the manuscript

`paper.tex` is a standalone source requiring a normal LaTeX distribution with
the packages named in its preamble. Run `pdflatex -interaction=nonstopmode
-halt-on-error paper.tex` twice in an output directory. The accompanying
`paper.md` is the accessible scientific manuscript; Pandoc was used to convert
its mathematical content to the supplied standalone source. LaTeX compilation
is a presentation step, not a scientific certificate replay.

The committed PDF is visually inspected page by page and checked for raw TeX.
Its exact bytes are bound by MANIFEST.sha256 and the publication asset ledger.
