import logging

def setlogging(config=None):
    """

    """
    if config == None:
        config = config_options() # sets default values
    levels = {
        'debug': DEBUG,
       'info': INFO
        }

    level = levels.get(config['log_level'])
    log_format = config['log_format']
    datefmt = config['log_datefmt']

    basicConfig(
        level   = level,
        format  = log_format,
        datefmt = datefmt)
