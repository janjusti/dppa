### Deep Protein Polarity Analyser (DPPA)
### Copyright (c) 2019 Jan Marans Agnella Justi & Mariana de França Costa

from .src import _core_methods

_core_methods_handler = _core_methods.CoreMethods()

def _main():
    _core_methods_handler.start_main()

def run(target_fn):
    return _core_methods_handler.start_run(target_fn)

def export(report_name, report_type, results_df_list):
    _core_methods_handler.start_export(report_name, report_type, results_df_list)

def set_debug_mode(isActive):
    _core_methods_handler.set_debug_mode(isActive)