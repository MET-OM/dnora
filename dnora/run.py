from subprocess import Popen

# =============================================================================
# RUN THE WAVE MODEL
# =============================================================================


def run_SWAN(input_file, swan_directory):
    print('Running SWAN----------------------->>>>>>>>>>>>>>>>>>>>>>>>>>')
    p = Popen(['./swanrun', '-input', input_file], cwd=swan_directory)
    p.wait()
