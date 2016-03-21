import logging.config


def config(debug=True):
    _ansi = {
        'reset': '\x1b[0m',
        'norm': '\x1b[22m',
        'purple': '\x1b[35m',
        'white': '\x1b[37m',
        'cyan': '\x1b[36m',
        'red': '\x1b[31;1m',
    }

    _1 = _ansi['purple']
    _2 = _ansi['purple']
    _h = _ansi['cyan']
    _t = _ansi['white']
    log_level = 'DEBUG' if debug else 'INFO'
    config_dict = {
        'version': 1,
        'disable_existing_loggers': False,
        'formatters': {
            'custom': {
                'format': (_1 + '[' + _2 + '%(asctime)s' +
                           _1 + ' - ' + _h + '%(name)s' +
                           _1 + ':' + _h + '%(lineno)s' +
                           _1 + ' - ' + _2 + '%(levelname)s' +
                           _1 + '] ' + _t + '%(message)s' +
                           _1 + _ansi['reset']
                           ),
                'datefmt': '%y-%m-%d %H:%M:%S',
            },
        },
        'filters': {
        },
        'handlers': {
            'null': {
                'level': 'DEBUG',
                'class': 'logging.NullHandler',
            },
            'console': {
                'level': 'DEBUG',
                'class': 'logging.StreamHandler',
                'formatter': 'custom'
            },
        },
        'loggers': {
            'root': {
                'handlers': ['null'],
                'propogate': True,
                'level': log_level,
            },
            '': {
                'handlers': ['console'],
                'propogate': True,
                'level': log_level,
            }
        }
    }
    logging.config.dictConfig(config_dict)
