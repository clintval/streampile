# Developer Documentation

## Local Setup

Install [uv](https://docs.astral.sh/uv/), then install the library and its development tools with:

```console
uv sync --locked
```

## Local Testing

To ensure all tests pass, run:

```console
uv run poe check-tests
```

The tests run the examples in the README too, so keep them working.

## Local Linting and Formatting

To check the lockfile, project metadata, format, lint, and types of all the code, and run every test, run:

```console
uv run poe check-all
```

To fix what can be fixed automatically, run:

```console
uv run poe fix-all
```

## Locking

`[tool.uv]` in `pyproject.toml` ignores releases younger than a week, and `uv.lock` records that setting.
Lock with the uv version pinned as `UV_VERSION` in [`tests.yml`](.github/workflows/tests.yml), and keep any user-level uv configuration out of the lock, or CI will find the lockfile stale:

```console
XDG_CONFIG_HOME="$(mktemp -d)" uvx uv@0.12.20 lock
```

## Commits

Commit titles follow [Conventional Commits](https://www.conventionalcommits.org), which group the release notes that [git-cliff](https://git-cliff.org) generates from `pyproject.toml`.
