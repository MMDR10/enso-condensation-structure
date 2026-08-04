import markdown, pathlib, subprocess

def md2pdf(md_path, pdf_path):
    md_text = pathlib.Path(md_path).read_text()
    html_body = markdown.markdown(md_text, extensions=['tables', 'fenced_code'])
    html = """<!DOCTYPE html>
<html><head><meta charset="utf-8"><style>
body { font-family: 'Georgia','Times New Roman',serif; max-width: 800px; margin: 40px auto; padding: 0 20px; line-height: 1.55; color: #111; font-size: 12pt; }
h1 { font-size: 17pt; line-height: 1.3; }
h2 { font-size: 14pt; margin-top: 1.6em; border-bottom: 1px solid #999; padding-bottom: 4px; }
h3 { font-size: 12.5pt; }
table { border-collapse: collapse; width: 100%; font-size: 9.5pt; margin: 1em 0; }
th, td { border: 1px solid #555; padding: 5px 8px; text-align: left; vertical-align: top; }
th { background: #eee; }
code { font-family: 'Courier New',monospace; font-size: 10pt; background:#f5f5f5; }
pre { background:#f5f5f5; padding:8px; }
blockquote { border-left: 3px solid #999; margin-left: 0; padding-left: 1em; color:#333; }
@page { margin: 2.2cm; }
</style></head><body>""" + html_body + "</body></html>"
    html_path = str(pdf_path).replace('.pdf', '.html')
    pathlib.Path(html_path).write_text(html)
    r = subprocess.run(['chromium', '--headless', '--disable-gpu', '--no-sandbox',
                        '--print-to-pdf=' + str(pdf_path), html_path],
                       capture_output=True, text=True)
    print(pdf_path, '->', r.returncode)

md2pdf('paper/paper_EN.md', 'paper/paper_EN.pdf')
md2pdf('paper/paper_ZH.md', 'paper/paper_ZH.pdf')
