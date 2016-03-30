'''
Created on Feb 20, 2015

@author: biiremployee

NOTE: _ANSI color coding will break if the color happens to be at a line break
'''

_ANSI = {'reset': '\x1b[0m',
         'norm': '\x1b[22m',
         'purple': '\x1b[35m',
         'white': '\x1b[37m',
         'cyan': '\x1b[36m',
         'red': '\x1b[31;1m',
         }

_NEWLINE = _ANSI['reset'] + "\n"
_PRIMARY = _ANSI['purple']
_SECONDARY = _ANSI['purple']
_HIGHLIGHT = _ANSI['cyan']
_TEXT = _ANSI['white']

_DEBUG = True
_LOG_LEVEL = 'DEBUG' if _DEBUG else 'INFO'
_LOGGING = {
    'version': 1,
    'disable_existing_loggers': False,
    'formatters': {
        'verbose': {
            'format': '%(levelname)s %(asctime)s %(module)s %(process)d ' +
                      '%(thread)d %(message)s'
        },
        'custom': {
            'format': (_PRIMARY + '[' +
                       _SECONDARY + '%(asctime)s' +
                       _PRIMARY + ' - ' +
                       _HIGHLIGHT + '%(name)s' +
                       _PRIMARY + ' - ' +
                       _HIGHLIGHT + '%(lineno)s CL' +
                       _PRIMARY + ' - ' +
                       _SECONDARY + '%(levelname)s' +
                       _PRIMARY + '] __beg__' +
                       _NEWLINE +
                       _TEXT + '%(message)s' +
                       _NEWLINE +
                       _PRIMARY + '__end__' +
                       _ANSI['reset']
                       )
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
        'django': {
            'handlers': ['null'],
            'propagate': True,
            'level': 'INFO',
        },
        'root': {
            'handlers': ['null'],
            'propogate': True,
            'level': _LOG_LEVEL,
        },
        'rempipe': {
            'handlers': ['null'],
            'propogate': True,
            'level': _LOG_LEVEL,
        },
        '': {
            'handlers': ['console'],
            'propogate': True,
            'level': _LOG_LEVEL,
        }
    }
}

import logging
import logging.config

_LOGGER_CONFIG_DONE = False


def config():
    global _LOGGER_CONFIG_DONE
    global _LOGGING
    if not _LOGGER_CONFIG_DONE:
        print('WARNING! Initializing custom logging')
        logging.config.dictConfig(_LOGGING)
        _LOGGER_CONFIG_DONE = True
