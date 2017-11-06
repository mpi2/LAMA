import nose
from nose import with_setup
from stats import automated_annotation
import numpy as np
import tempfile
from os.path import join, abspath, dirname
import SimpleITK as sitk
import numpy as np


def make_test_volumes():
	"""
	Make a test volume. 9 labels
	"""
	labels = np.zeros(shape=(90,90,90))
	for i, x in enumerate(range(0, 90, 10)):
		labels[x: x+90, x:x+90, x:x+90] = i

	np.random.seed(44)

	label_info ="""################################################
# ITK-SnAP Label Description File
# File format:
# IDX   -R-  -G-  -B-  -A--  VIS MSH  LABEL
# Fields:
#    IDX:   Zero-based index
#    -R-:   Red color component (0..255)
#    -G-:   Green color component (0..255)
#    -B-:   Blue color component (0..255)
#    -A-:   Label transparency (0.00 .. 1.00)
#    VIS:   Label visibility (0 or 1)
#    IDX:   Label mesh visibility (0 or 1)
#  LABEL:   Label description
#  TERM:    Ontology term eg. emapa # added by Neil
################################################
 0     0    0    0      0  0  0    'Clear Label'		        emapa1
1   223  169   45       1  0  0    'brain_cerebral_aqueduct'	emapa2
2   194    4  114       1  0  0    'brain_fourth_ventricle'	    emapa3
3   246  178  234       1  0  0    'brain_hypothalamus_left'	emapa4
4   95  143  242        1  0  0    'brain_hypothalamus_right'	emapa5
5   10  169   45        1  0  0    'brain_cerebral_thing'	    emapa6
6   20    4  114        1  0  0    'brain_fourth_thing'		    emapa7
7   30  178  234        1  0  0    'brain_hypothalamus_thing'	emapa8
8   40  143  242        1  0  0    'brain_hypothalamus_thing'	emapa9
"""

	# Write out the labelmap
	labelmap_file = tempfile.NamedTemporaryFile('r+', suffix='.nrrd')
	label_img = sitk.GetImageFromArray(labels)
	sitk.WriteImage(label_img, labelmap_file.name)


	label_info_file = tempfile.NamedTemporaryFile()
	outfile = join(dirname(abspath(__file__)), 'test_output', 'autoannotator_result.csv')





	# Write out the label info file
	with open(label_info_file.name, 'w') as lif:
		lif.write(label_info)

	# Here I'm using the labelmap as the dummy t-stats file as well. It should produce top hits for the
	#  higher -labelled organs
	ann = automated_annotation.Annotator(labelmap_file.name, label_info_file.name, labelmap_file.name, outfile)
	ann.annotate()

def test_with_volume():
	pass