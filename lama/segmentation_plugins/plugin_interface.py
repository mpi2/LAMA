from logzero import logger as logging
from lama.registration_pipeline.validate_config import LamaConfig
from lama.common import cfg_load, write_array
import importlib.util
from lama.elastix import INVERT_CONFIG


def secondary_segmentation(config: LamaConfig):
    """
    Use user-added scripts to segment/cleanup organs

    Parameters
    ----------
    config

    Returns
    -------

    """

    plugin_dir = config.config_dir / config['seg_plugin_dir']

    if not plugin_dir.is_dir():
        logging.error(f'Cannot find plugin director: {plugin_dir}')
        return

    # Find the directories containing the segmentations
    # Get the final inversion stage
    invert_config = config['inverted_transforms'] / INVERT_CONFIG
    segmentation_dir = cfg_load(invert_config)['inversion_order'][-1] # rename to segmentation stage
    inverted_label_dir = config['inverted_labels'] / segmentation_dir
    initial_segmentation_path = next(inverted_label_dir.glob('**/*.nrrd'))

    first_reg_dir = config['root_reg_dir'] / config['registration_stage_params'][0]['stage_id']  # usually rigid
    image_to_segment = next(first_reg_dir.glob('**/*.nrrd'))

    segmentations = []

    for plugin_src in [x for x in plugin_dir.iterdir() if str(x).endswith('.py') and x.name != 'plugin_interface.py']:

        # catch all exceptions as we don't want plugin crashing the pipeline
        try:
            spec = importlib.util.spec_from_file_location(plugin_src.stem,  str(plugin_src))
            plugin = importlib.util.module_from_spec(spec)
            spec.loader.exec_module(plugin)

            new_segmetation = plugin.run(image_to_segment, initial_segmentation_path)

        except Exception as e:
            logging.error(f'Plugin {plugin_src} failed\n{e}')
        else:
            segmentations.append(new_segmetation)

    if not segmentations:
        logging.error(f'No segmentations returned from {plugin_src.name}')

    # Merge all the segmentations into a single label map. If there are any overlaps, the plugin called last will have
    # priority

    seg = None

    for s in segmentations:
        if not seg:
            seg = s
            continue
        seg[s != 0] = s[s != 0]

    additional_seg_dir = config.mkdir('additional_seg_dir')
    write_array(seg, additional_seg_dir / f'{config.config_dir.name}_additonal_seg.nrrd')   # TODO include specimen name
