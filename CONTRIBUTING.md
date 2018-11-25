# Contributing to proEQUIB

The following guidelines are designed for contributors to the proEQUIB IDL package, which is 
hosted in the [proEQUIB repository](https://github.com/equib/proEQUIB) on GitHub. 

## Reporting Issues

The [issue tracker](https://github.com/equib/proEQUIB/issues) is used to report bugs, request new functionality, and discuss improvements. 
For reporting a bug or a failed function or requesting a new feature, you can simply open an issue 
in the [issue tracker](https://github.com/equib/proEQUIB/issues) of the 
[proEQUIB repository](https://github.com/equib/proEQUIB). If you are reporting a bug, please also include a minimal code
example that reproduces the problem, and IDL version you are using.

## Contributing Code

Fo contributing code to proEQUIB, you need to set up your [GitHub](https://github.com) 
account if you do not have and sign in, and request your change(s) or contribution via 
opening a pull request against the ``master``
branch in your fork of the [proEQUIB repository](https://github.com/equib/proEQUIB). 

To contribute to this package, you need to follow these steps:

- Open a new issue for new feature or failed function in the [Issue tracker](https://github.com/equib/proEQUIB/issues).
- Fork the [proEQUIB repository](https://github.com/equib/proEQUIB) to your GitHub account.
- Clone your fork of the [proEQUIB repository](https://github.com/equib/proEQUIB):

      $ git clone git@github.com:your-username/proEQUIB.git
      
- Make your change(s) in the `master` branch of your cloned fork.
- Make sure that it passes all tests and there is no error.
- Push yout change(s) to your fork in your GitHub account.
- [Submit a pull request][pr], mentioning what issue has been addressed.

[pr]: https://github.com/equib/proEQUIB/compare/

Then, you are waiting, until your contribution is checked and merged into the original repository. 
We will contact you if there is a problem in your code.

While you are opening a pull request for your contribution, be sure that you have included:

* **Code** which you are contributing to this package.

* **Documentation** of this code if it provides new functionality. This should be a
  description of new functionality added to the API documentation (in ``docs/``). 

- **Tests** of this code to make sure that the previously failed function or the new functionality now works properly.

- **Revision history** if you fixed a bug in the previously failed function or add a code for new functionality, you should
well document your change(s) or addition in the *Revision History* entry of the changed or added function in your code.
