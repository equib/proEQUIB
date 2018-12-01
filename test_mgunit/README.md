## Test with mgunit

To test this package with [mgunit](https://github.com/mgalloy/mgunit), you need to install ``mgunit`` on your machine:

    git clone --recursive https://github.com/mgalloy/mgunit.git

The test files need to have accessible to AtomNeb database in ``externals/atomneb/atomic-data`` and ``externals/atomneb/atomic-data-rc``. To have all of them in ``externals`` directory of proEQUIB, you can download proEQUIB with all dependencies as follows 

    git clone --recursive https://github.com/equib/proEQUIB.git

To run the test, you need to run the following command:

    idl test_all.pro

This will produce two files: ``test-results.log`` and ``test-results.html``. All the functions should pass the test without any failure.



