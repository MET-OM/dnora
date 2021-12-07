# =============================================================================
# IMPORT dnora
# =============================================================================
import sys
dnora_directory = '../'
sys.path.insert(0, dnora_directory)
from dnora import trg

trgrid = trg.TrGrid(name='test')

trgrid.import_triang(trg.read_tr.SmsReader(filename='WW3_grid_01.2dm'))
