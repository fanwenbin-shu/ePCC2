import os
import time
import logging

def check_os():
    '''Check operating system (OS) for a better performance. '''
    os_name = os.name.lower()
    if os_name != 'posix':
        if os_name == 'nt':
            print('Do not support Windows OS. ', end='')
        print('Support LINUX only! ')
        exit()
    return

def print_header():

    logging.info('Polarized Crystal Charge (PCC)')
    logging.info('Version : 2.0')
    logging.info('')
    logging.info('Program Author : Wenbin FAN (fanwenbin@shu.edu.cn)')
    logging.info('                 Feng XU (xufeng13@shu.edu.cn)')
    logging.info('Supervisor : Yongle LI (yongleli@shu.edu.cn)')
    logging.info('Thank Xiaoqing ZHU for the help of mastering PCC! ')
    logging.info('')

    return

def log_format():
    try:
        import logging
    except ImportError:
        raise ImportError('Package `logging` is necessary. Please Check it. ')

    logFormatter = logging.Formatter("%(asctime)s [%(levelname)-5.5s] %(message)s")
    logger = logging.getLogger()
    logger.setLevel(logging.DEBUG)

    fileHandler = logging.FileHandler("PyPCC.log")
    fileHandler.setFormatter(logFormatter)
    logger.addHandler(fileHandler)

    consoleHandler = logging.StreamHandler()
    consoleHandler.setFormatter(logging.Formatter("%(message)s"))
    logger.addHandler(consoleHandler)
    return

def divmod_excel(n):
    a, b = divmod(n, 26)
    if b == 0:
        return a - 1, b + 26
    return a, b

def get_alphabet(num):
    chars = []
    num += 1
    while num > 0:
        num, d = divmod_excel(num)
        chars.append(chr(d - 1 + 65))
    return ''.join(reversed(chars))
