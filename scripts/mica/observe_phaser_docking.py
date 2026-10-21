import faulthandler, runpy, sys, time
faulthandler.enable()
faulthandler.dump_traceback_later(300, repeat=True)
import phenix.autosol.run_phaser_mr as mr
original = mr.runMR_AUTO
def observed(inp):
    print("DIAGNOSTIC: entering runMR_AUTO", time.time(), flush=True)
    inp.setMUTE(False)  # Output verbosity only; retain every search setting.
    result = original(inp)
    print("DIAGNOSTIC: leaving runMR_AUTO", time.time(), flush=True)
    return result
mr.runMR_AUTO = observed
sys.argv = sys.argv[1:]
runpy.run_path(sys.argv[0], run_name='__main__')
