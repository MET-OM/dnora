from subprocess import Popen

# =============================================================================
# RUN THE WAVE MODEL
# =============================================================================


def run_SWAN(input_file, swan_directory):
    print('Running SWAN----------------------->>>>>>>>>>>>>>>>>>>>>>>>>>')
    p = Popen(['./swanrun', '-input', input_file], cwd=swan_directory)
    p.wait()


def run_SWASH(input_file, swash_directory):
    print('Running SWASH----------------------->>>>>>>>>>>>>>>>>>>>>>>>>>')
    p = Popen(['./swashrun', '-input', input_file], cwd=swash_directory)
    p.wait()
