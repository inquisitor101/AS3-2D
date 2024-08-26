# Contributing Guidelines

Collaboration is a cornerstone of this project. We warmly welcome and highly encourage every contribution, no matter the size. Your involvement is deeply appreciated and vital to the success of the project.

## Contents

* [Code of Conduct](#book-code-of-conduct)
* [Discussions](#speech_balloon-discussions)
* [Coding Style](#nail_care-coding-style)
* [Opening an Issue](#monocle_face-opening-an-issue)
* [Request a New Feature](#bulb-request-a-new-feature)
* [Submit a Pull Request](#pleading_face-submit-a-pull-request)
* [Code Review](#person_fencing-code-review)


# :book: Code of Conduct

Please review and adhere to the [Code of Conduct](./CODE_OF_CONDUCT.md) at all times. It is essential for maintaining a respectful and collaborative environment.


# :speech_balloon: Discussions

Open and constructive discussions are key to a successful and collaborative project, especially for open-source initiatives like AS3-2D. We highly encourage you to share ideas, provide feedback and engage in discussions through [GitHub Discussions](https://docs.github.com/en/discussions/quickstart). Additionally, you can join our public [Slack channel](https://join.slack.com/t/as3-2d/shared_invite/zt-2nxm0hq2u-vwV9I8wIru1YlkN9sUnhQA) for more informal conversations.


# :nail_care: Coding Style

Maintaining a consistent coding style is crucial for keeping the project organized and manageable. Please follow the [AS3-2D coding guidelines](./coding_style.md) to ensure uniformity across the codebase.


# :monocle_face: Opening an Issue

Reporting issues or bugs is greatly appreciated and provides a valuable opportunity to contribute to the project. We use [GitHub Issues](https://docs.github.com/en/issues/tracking-your-work-with-issues/about-issues) for bug management, which you can learn about from [here](https://help.github.com/en/github/managing-your-work-on-github/creating-an-issue). Please keep the following guidelines in mind:

* Before reporting a new bug, please review the existing issues to avoid creating duplicates. 
* Ensure you fill out the bug report [template](./ISSUE_TEMPLATE/bug_report.md) completely.


# :bulb: Request a New Feature

Adding new features is encouraged and always welcome. However, please keep in mind that some features may be beyond the current scope or timeline of AS3-2D. That said, we will carefully review your proposal. A few important points to consider:

* Avoid submitting duplicate feature requests.
* Ensure that the feature request [template](./ISSUE_TEMPLATE/feature_request.md) is fully completed.


# :pleading_face: Submit a Pull Request

If you're reading this, you're likely looking for guidance on how to submit a pull request (PR). We use the [fork-and-pull](https://github.com/susam/gitpr) workflow, which is outlined below. Before diving into your PR, we recommend opening an issue or starting a discussion about your idea. This allows you to gather feedback and possibly find collaborators to help solve the problem.

## To start contributing to AS3-2D, follow these steps:
1. Fork the repository to your own GitHub account.
2. Clone your forked repository to your local machine.
3. Set up the upstream repository.
4. Create a new branch for your changes.
5. Commit your changes locally and push them to your forked repository.
6. Once your feature is complete, open a PR using the [template](./pull_request_template.md).
7. If your PR needs adjustments, repeat steps 5 and 6 until it is approved.
8. (Optional) Delete your branch after your PR is approved.

## Example Walkthrough: Contributing to the Repository

Follow these steps to contribute to the project:

### Step 1: Fork the Original Repository

1. Open the original repository (referred to as **upstream**) in your browser.
2. Click the "Fork" button to create a copy of the repository under your GitHub account.

### Step 2: Clone Your Forked Repository

```bash
git clone https://github.com/YourGithubAccount/AS3-2D.git
```
This command downloads your forked repository to your local machine.

### Step 3: Add the Original Repository as a Remote

```bash
git remote add upstream https://github.com/inquisitor101/AS3-2D.git
git remote -v
```
This step links your local repository to the original repository, allowing you to pull in updates.

### Step 4: Create a New Feature Branch
```bash
git checkout -b feature_name
```
Create and switch to a new branch where you’ll make your changes.

### Step 5: Make Your Changes and Push the Branch
```bash
git add .
git commit -m "Description of the changes made."
git push origin feature_name
```
Stage your changes, commit them with a descriptive message, and push the branch to your forked repository on GitHub.

### Step 6: Open a Pull Request
1. Go to your forked repository on GitHub using your browser.
2. Click the "Compare & pull request" button.
3. Follow the template instructions when submitting your pull request (PR).

By following these steps, you'll contribute your changes to the project effectively.


# :person_fencing: Code Review

* Focus on the code, not the author: remember, the goal is to review and improve the code, not to critique the individual who wrote it.
* Constructive criticism: whether giving or receiving feedback, aim to be constructive and helpful. This approach fosters growth and improvement for everyone involved.
* Be respectful: mistakes are a natural part of the process. Approach reviews with kindness and understanding.
* Follow guidelines: please adhere to the guidelines outlined in this document to ensure a consistent and effective review process.



