# Generic tests file for use with the pytest module

# Test functions have names starting with "test_"
# These will all be discovered by simply running pytest from the repo root.

# Tests can run python code for the package, or invoke command-line usage.

# Writing good tests is as much an art as a science.
# As the software gets more complex, more complex testing methods are needed.

# But generally, the best tests follow this structure:
# Arrange - set up conditions for the test
# Act - call some function, method or command
# Assert - check that some invariant is true.

import os
    
# This simple, system-level test just runs the command line version of the package as in the example script.
# "monkeypatch" is a text fixture that lets you change context scoped to the test.
def test_system_cli(monkeypatch):
    monkeypatch.chdir("./examples/scripts")
    exitcode = os.system("./run_example.sh")
    assert exitcode == 0, "exitcode != 0"