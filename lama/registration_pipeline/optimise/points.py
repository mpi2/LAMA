"""
Convert between points fieltypes needed for the transforming of points.

Examples
--------
# Add points file frrom elastix
p = Points().set_elastix_output('output_points.txt')

# Save as Slicer format
p.write_slicer(Path('points.fcsv'))

# Get points in a pandas DataFrame with columns x,y,z
df = p.points

# Or do all the above in one line
df = Points().set_elastix_output('output_points.txt').write_slicer('points.fcsv').points
"""


import pandas as pd
import re


class Points:

    def __init__(self):
        self.points = None

    def set_slicer(self, points_file):
        slicer_col_row_prefix = '# columns = '
        with open(points_file, 'r') as fh:
            for line in fh:
                if line.startswith(slicer_col_row_prefix):
                    header = line.strip(slicer_col_row_prefix).strip().split(',')
                    break

        df = pd.read_csv(points_file, comment='#', header=None)
        df.columns = header
        df = df[['x', 'y', 'z']]
        self.points = df

        return self

    def set_elastix_output(self, points_file):
        result = []
        with open(points_file, 'r') as fh:
            for line in fh:
                index_str = line.split(';')[3]
                arraystr = re.search('\[(.+)\]', index_str)
                if arraystr:
                    a = [float(x) for x in arraystr.group(1).strip().split(' ')]
                    result.append(a)
        df = pd.DataFrame.from_records(result)
        df.columns = ['x', 'y', 'z']
        self.points = df

        return self

    def write_elastix_input(self, output_file: str):

        if self.points is None:
            raise ValueError('No points have been set')

        with open(output_file, 'w') as fh1:
            # Write the header
            fh1.write(f'index\n{len(self.points)}\n')

            # write the points
            for i, row in self.points.iterrows():
                fh1.write(f'{row.x} {row.y} {row.z}\n')

        return self

    def write_slicer(self, output_file: str):
        """

        Parameters
        ----------
        output_file


        Notes
        -----

        Example output

        # Markups fiducial file version = 4.8
        # CoordinateSystem = 0
        # columns = id,x,y,z,ow,ox,oy,oz,vis,sel,lock,label,desc,associatedNodeID
        vtkMRMLMarkupsFiducialNode_0,93.2369,371.84,595.904,0,0,0,1,1,1,0,F-1,,vtkMRMLScalarVolumeNode1
        vtkMRMLMarkupsFiducialNode_1,304.938,371.84,588.265,0,0,0,1,1,1,0,F-2,,vtkMRMLScalarVolumeNode1

        """

        header = ("# Markups fiducial file version = 4.8\n"
        "# CoordinateSystem = 0\n"
        "# columns = id,x,y,z,ow,ox,oy,oz,vis,sel,lock,label,desc,associatedNodeID\n")

        row_template = 'vtkMRMLMarkupsFiducialNode_{},{},{},{},0,0,0,1,1,1,0,F-{},,vtkMRMLScalarVolumeNode{}\n'

        if self.points is None:
            raise ValueError('No points have been set')

        with open(output_file, 'w') as fh:
            # Bodge to get rid of negative values in the slicer file
            df = self.points.copy()
            df = df[['x', 'y', 'z']].abs()

            fh.write(header)
            for i, row in df.iterrows():
                fh.write(row_template.format(i, row.x, row.y, row.z, i, i))

        return self
