"""
Create grids of images for qc.
"""

from pathlib import Path
import shutil


class HtmlGrid:
    """
    Create a grid of images
    """

    def __init__(self, title):
        self.title = title
        padding = 20
        margin_bottom = 20

        self.html = ('<html>\n'
                     '<head>\n'
                     '<script src="https://code.jquery.com/jquery-3.5.1.js"></script>\n'  # For tooltips
                    # '<script src="https://code.jquery.com/ui/1.12.1/jquery-ui.js"></script>\n'
                     '</head>\n'
                     '<body>\n'
                     '<div class="title">{}</div>\n').format(title)
        self.num_rows = 0

    def next_row(self, title: str = ''):
        # End the old row
        if self.num_rows > 0:
            self.html += '</div>\n'

        # open a new row div
        self.html += '<div class="clear"></div>'
        self.html += f'<div class="row_title">{title}</div>\n'
        self.html += '<div class="row clear">\n'

    def next_image(self, img_path: Path, caption: str = '', tooltip: str = ''):
        self.html += ('<div data-tooltip class="image" title={}>\n'
                      '<img src="{}">\n'
                      '<div class="top-left">{}</div>'
                      '</div>').format(tooltip, img_path, str(caption))

    def save(self, outpath: Path):
        # close the final row
        self.html += '</div>\n'
        # close the file
        self.html += '</body></html>\n'
        with open(outpath, 'w') as fh:
            fh.write(self.html)
            fh.write(self.get_css())

    def get_css(self):
        css = ('<style>'
               'body{font-family: Arial}\n'
               '.title{width: 100%; padding: 2px; background-color: lightblue; margin-bottom: 5px}\n'
               '.row{auto;white-space: nowrap;}\n'
               '.image{display: inline-block;text-align: center;color: white; position:relative}\n'
               'img{white-space: nowrap;}'
               '.clear{clear: both}\n'
               '.row_label{font-size: 2em}\n'
               '.top-left{position: absolute;top: 8px;left: 16px;}'
               '</style>')
        return css

