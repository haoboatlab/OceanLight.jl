# Contributors Guide

Thank you for considering contributions to OceanLight! We hope this guide helps.

Feel free to ask us questions and chat with us at any time about any topic at all
by:

* [Opening a GitHub issue](https://github.com/haoboatlab/OceanLight.jl/issues/new)

## Creating issues

The simplest way to contribute to OceanLight is to create or comment on issues.

The most useful bug reports:

* Provide an explicit code snippet --- not just a link --- that reproduces the bug in the latest tagged version of OceanLight. This is sometimes called the ["minimal working example"](https://en.wikipedia.org/wiki/Minimal_working_example). Reducing bug-producing code to a minimal example can dramatically decrease the time it takes to resolve an issue.

* Paste the _entire_ error received when running the code snippet, even if it's unbelievably long.

* Use triple backticks (e.g., ````` ```some_code; and_some_more_code;``` `````) to enclose code snippets, and other [markdown formatting syntax](https://docs.github.com/en/github/writing-on-github/getting-started-with-writing-and-formatting-on-github/basic-writing-and-formatting-syntax) to make your issue easy and quick to read.

* Report the OceanLight version, Julia version, machine (especially if using a GPU) and any other possibly useful details of the computational environment in which the bug was created.

## But I want to _code_!

* New users help write OceanLight code and documentation by [forking the OceanLight repository](https://docs.github.com/en/github/collaborating-with-pull-requests/working-with-forks), [using git](https://guides.github.com/introduction/git-handbook/) to edit code and docs, and then creating a [pull request](https://docs.github.com/en/github/collaborating-with-pull-requests/proposing-changes-to-your-work-with-pull-requests/creating-a-pull-request-from-a-fork). Pull requests are reviewed by OceanLight collaborators.

* A pull request can be merged once it is reviewed and approved by collaborators. If the pull request author has write access, they have the responsibility of merging their pull request. Otherwise, OceanLight.jl collaborators will execute the merge with permission from the pull request author.

* Note: for small or minor changes (such as fixing a typo in documentation), the [GitHub editor](https://docs.github.com/en/github/managing-files-in-a-repository/managing-files-on-github/editing-files-in-your-repository) is super useful for forking and opening a pull request with a single click.

* Write your code with love and care. In particular, conform to existing OceanLight style and formatting conventions. For example, we love verbose and explicit variable names, use `TitleCase` for types, `snake_case` for objects, and always.put.spaces.after.commas. For formatting decisions we loosely follow the [YASGuide](https://github.com/jrevels/YASGuide). It's worth few extra minutes of our time to leave future generations with well-written, readable code.

## What is a "collaborator" and how can I become one?

* Collaborators have permissions to review pull requests and status allows a contributor to review pull requests in addition to opening them. Collaborators can also create branches in the main OceanLight repository.

* We ask that new contributors try their hand at forking OceanLight, and opening and merging a pull request before requesting collaborator status.

## What's a good way to start developing OceanLight?

* Try to run OceanLight and play around with it. If you run into any problems or find it difficult
  to use or understand, please open an issue!

* Write up an example or tutorial on how to do something useful with
  OceanLight, like how to set up a new configuration.

* Improve documentation or comments if you found something hard to use.

* Implement a new feature if you need it to use OceanLight.

If you're interested in working on something, let us know by commenting on existing issues or 
by opening a new issue. This is to make sure no one else is working on the same issue and so 
we can help and guide you in case there is anything you need to know beforehand.

## Ground Rules

* Each pull request should consist of a logical collection of changes. You can
  include multiple bug fixes in a single pull request, but they should be related.
  For unrelated changes, please submit multiple pull requests.

* Do not commit changes to files that are irrelevant to your feature or bugfix
  (eg: `.gitignore`).

* Be willing to accept criticism and work on improving your code; we don't want
  to break other users' code, so care must be taken not to introduce bugs. We
  discuss pull requests and keep working on them until we believe we've done a
  good job.

* Be aware that the pull request review process is not immediate, and is
  generally proportional to the size of the pull request.

## Reporting a bug

The easiest way to get involved is to report issues you encounter when using
OceanLight or by requesting something you think is missing.

* Head over to the [issues](hhttps://github.com/haoboatlab/OceanLight.jl/issues) page.

* Search to see if your issue already exists or has even been solved previously.

* If you indeed have a new issue or request, click the "New Issue" button.

* Please be as specific as possible. Include the version of the code you were using, as
  well as what operating system you are running. The output of Julia's `versioninfo()`
  and `] status` is helpful to include. Try your best to include a complete, ["minimal working example"](https://en.wikipedia.org/wiki/Minimal_working_example) that reproduces the issue.

## Setting up your development environment

* Install [Julia](https://julialang.org/) on your system.

* Install `git` on your system if it is not already there (install XCode command line tools on
  a Mac or `git bash` on Windows).

* Login to your GitHub account and make a fork of the
  [OceanLight repository](https://github.com/haoboatlab/OceanLight.jl) by
  clicking the "Fork" button.

* Clone your fork of the OceanLight repository (in terminal on Mac/Linux or git shell/
  GUI on Windows) in the location you'd like to keep it.
  ```
  git clone https://github.com/haoboatlab/OceanLight.jl.git
  ```

* Navigate to that folder in the terminal or in Anaconda Prompt if you're on Windows.

* Connect your repository to the upstream (main project).
  ```
  git remote add OceanLight https://github.com/haoboatlab/OceanLight.jl.git
  ```

* Create the development environment by opening Julia via `julia --project` then
  typing in `] instantiate`. This will install all the dependencies in the Project.toml
  file.

* You can test to make sure OceanLight works by typing in `] test`. Doing so will run all
  the tests (and this can take a while).

Your development environment is now ready!

## Pull Requests

We follow the [ColPrac guide](https://github.com/SciML/ColPrac) for collaborative practices.
We ask that new contributors read that guide before submitting a pull request.

Changes and contributions should be made via GitHub pull requests against the ``main`` branch.

When you're done making changes, commit the changes you made. Chris Beams has written a 
[guide](https://chris.beams.io/posts/git-commit/) on how to write good commit messages.

When you think your changes are ready to be merged into the main repository, push to your fork
and [submit a pull request](https://github.com/haoboatlab/OceanLight.jl/compare/).

**Working on your first Pull Request?** You can learn how from this _free_ video series
[How to Contribute to an Open Source Project on GitHub](https://egghead.io/courses/how-to-contribute-to-an-open-source-project-on-github), Aaron Meurer's [tutorial on the git workflow](https://www.asmeurer.com/git-workflow/), or the guide [“How to Contribute to Open Source"](https://opensource.guide/how-to-contribute/).

## Documentation

Now that you've made your awesome contribution, it's time to tell the world how to use it.
Writing documentation strings is really important to make sure others use your functionality
properly. Didn't write new functions? That's fine, but be sure that the documentation for
the code you touched is still in great shape. It is not uncommon to find some strange wording
or clarification that you can take care of while you are here.

You can preview how the Documentation will look like after merging by building the documentation 
locally. From the main directory of your local repository call

```
JULIA_DEBUG=Documenter julia --project=docs/ docs/make.jl
```

and then open `docs/build/index.html` in your favorite browser. Providing the environment variable 
`JULIA_DEBUG=Documenter` will provide with more information in the documentation build process and
thus help figuring out a potential bug.

## Credits

This contributor's guide is heavily based on the excellent [OceanBioMe contributors guide](https://oceanbiome.github.io/OceanBioME.jl/stable/contributing/) and [Oceananigans contributors guide](https://clima.github.io/OceananigansDocumentation/stable/contributing/) which in turn is based on the [MetPy contributor's guide](https://github.com/Unidata/MetPy/blob/main/CONTRIBUTING.md).