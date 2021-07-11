from PCC.main import ePCC
from PCC.mist import log_format

if __name__ == '__main__':
    log_format()

    ins = ePCC()
    ins.main()