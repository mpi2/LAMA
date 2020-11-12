"""
Create grids of images for qc.
"""

from pathlib import Path
import shutil
import base64


class HtmlGrid:
    """
    Create a grid of images
    """

    def __init__(self, title, single_file=False):
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
        self.single_file = single_file

    def next_row(self, title: str = ''):
        # End the old row
        if self.num_rows > 0:
            self.html += '</div>\n'

        # open a new row div
        self.html += '<div class="clear"></div>'
        self.html += f'<div class="row_title">{title}</div>\n'
        self.html += '<div class="row clear">\n'

    def next_image(self, img_path: Path, caption: str = '', tooltip: str = ''):
        # data_uri = base64.b64encode(open('Graph.png', 'rb').read()).decode('utf-8')
        # img_tag = '<img src="data:image/png;base64,{0}">'.format(data_uri)
        # print(img_tag)
        if self.single_file:
            with open(img_path , 'rb') as fh:
                data_uri = base64.b64encode(fh.read()).decode('ascii')

            img_tag = '<img src="data:image/png;base64,{0}" >'.format(data_uri)
        else:
            img_tag = f'<img src="{str(img_path)}">'

        self.html += ('<div data-tooltip class="image" title={}>\n'
                      '{}\n'
                      '<div class="top-left">{}</div>'
                      '</div>').format(tooltip, img_tag, str(caption))

    def save(self, outpath: Path, height=600, width=None):
        # close the final row
        self.html += '</div>\n'
        # close the file
        self.html += '</body></html>\n'
        with open(outpath, 'w') as fh:
            fh.write(self.html)
            fh.write(self.get_css(height, width))

    def get_css(self, height=300, width=None):
        if width:
            how = f'width:{width}'
        else:
            how = f'height:{height}'


        css = ('<style>'
               'body{{font-family: Arial}}\n'
               '.title{{width: 100%; padding: 2px; background-color: lightblue; margin-bottom: 5px}}\n'
               '.row{{auto;white-space: nowrap;}}\n'
               '.image{{display: inline-block;text-align: center;color: white; position:relative}}\n'
               # 'img{white-space: nowrap; height: 100%; width: auto; max-width: MaxSize; max-height: MaxSize;}'
                'img{{white-space: nowrap; {}}}\n'
               '.clear{{clear: both}}\n'
               '.row_label{{font-size: 2em}}\n'
               '.top-left{{position: absolute;top: 8px;left: 16px;}}'
               '</style>').format(how)
        return css

