# Lint-only option from https://stackoverflow.com/a/60045845/12968623.
def pytest_collection_modifyitems(session, config, items):
    if config.getoption("--lint-only"):
        lint_items = []
        for linter in ["black", "flake8"]:
            if config.getoption(f"--{linter}"):
                lint_items.extend(
                    [item for item in items if item.get_closest_marker(linter)]
                )
        items[:] = lint_items


def pytest_addoption(parser):
    parser.addoption(
        "--lint-only",
        action="store_true",
        default=False,
        help="Run linting checks only.",
    )
