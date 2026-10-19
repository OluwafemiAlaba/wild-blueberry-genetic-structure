from pathlib import Path
import base64
from reportlab.pdfgen import canvas
from reportlab.lib.utils import ImageReader
from reportlab.pdfbase import pdfmetrics
from reportlab.pdfbase.ttfonts import TTFont
from matplotlib.font_manager import findfont
pdfmetrics.registerFont(TTFont('FigureArial',findfont('DejaVu Sans')))
from matplotlib.font_manager import FontProperties
pdfmetrics.registerFont(TTFont('FigureArialBold',findfont(FontProperties(family='DejaVu Sans',weight='bold'))))

ROOT = Path(__file__).resolve().parents[2]
OUT = Path(__file__).resolve().parent
SRC = ROOT / 'inputs/manuscript_media/image3.png'
W, H = 1320, 730
pdf = canvas.Canvas(str(OUT / 'Figure2_original_relabelled.pdf'), pagesize=(W,H))
pdf.setFillColorRGB(1,1,1); pdf.rect(0,0,W,H,fill=1,stroke=0)
svg = [f'<svg xmlns="http://www.w3.org/2000/svg" xmlns:xlink="http://www.w3.org/1999/xlink" width="{W}" height="{H}" viewBox="0 0 {W} {H}"><rect width="100%" height="100%" fill="white"/>']
encoded = base64.b64encode(SRC.read_bytes()).decode()

def original_clip(name, sx, sy, width, height, dx, dy):
    # Embed the original raster without changing its plotted pixels or colors.
    pdf.saveState(); p=pdf.beginPath(); p.rect(dx,H-dy-height,width,height)
    pdf.clipPath(p,stroke=0,fill=0)
    pdf.drawImage(ImageReader(str(SRC)),dx-sx,H-(dy-sy)-780,width=1140,height=780)
    pdf.restoreState()
    svg.append(f'<defs><clipPath id="{name}"><rect x="{dx}" y="{dy}" width="{width}" height="{height}"/></clipPath></defs>')
    svg.append(f'<image x="{dx-sx}" y="{dy-sy}" width="1140" height="780" clip-path="url(#{name})" xlink:href="data:image/png;base64,{encoded}"/>')

def text(x,y,s,size=25,anchor='middle',bold=False):
    pdf.setFillColorRGB(0,0,0); pdf.setFont('FigureArialBold' if bold else 'FigureArial',size)
    f={'middle':pdf.drawCentredString,'end':pdf.drawRightString,'start':pdf.drawString}[anchor]
    f(x,H-y,s)
    svg.append(f'<text x="{x}" y="{y}" font-family="Arial, Helvetica, sans-serif" font-size="{size}" text-anchor="{anchor}" font-weight="{"bold" if bold else "normal"}">{s}</text>')

original_clip('chromosomes',85,143,938,597,119,103)
for i in range(12): text(103,165+50*i-40+8,f'Chr{i+1}',25,'end')
pdf.setStrokeColorRGB(0,0,0);pdf.setLineWidth(2)
pdf.line(120,H-84,1052,H-84)
svg.append('<path d="M120 84H1052" stroke="black" stroke-width="2"/>')
ticks=[(86,'0'),(190,'6.4'),(293,'12.8'),(396,'19.2'),(500,'25.6'),(604,'32'),(707,'38.4'),(811,'44.8'),(914,'51.2'),(1018,'57')]
for x,label in ticks:
    x+=34; pdf.line(x,H-84,x,H-91)
    svg.append(f'<path d="M{x} 84v7" stroke="black" stroke-width="2"/>')
    text(x,68,label,23)
text(586,30,'Genomic position (Mb)',27)
original_clip('legend',1032,423,34,317,1120,285)
text(1172,242,'SNPs per',24,bold=True)
text(1172,270,'1-Mb window',24,bold=True)
for j,label in enumerate(['0','1','7','13','19','25','31','37','43','49','55']):
    text(1166,285+(j+.5)*(317/11)+8,label,25,'start')
svg.append('</svg>')
(OUT/'Figure2_original_relabelled.svg').write_text('\n'.join(svg))
pdf.showPage();pdf.save()
(OUT/'Figure2_original_relabelled_caption.txt').write_text('Figure 2. Distribution of SNP markers across the 12 chromosomes of the Vaccinium darrowii reference genome. Colors indicate the number of SNPs within a 1-Mb window; genomic positions are shown in megabases (Mb).\n')
