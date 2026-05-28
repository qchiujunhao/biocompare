# Contributing

`biocompare` is built around small comparator plugins that return a shared
`ConcordanceReport`. Contributions should keep that contract stable and add
tests that demonstrate the scientific behavior being compared.

## Development Setup

```bash
python3 -m unittest discover -s tests
```

The core package intentionally uses only the Python standard library. Optional
heavy dependencies should be isolated behind comparator-specific code paths or
third-party plugins.

## Comparator Checklist

- Implement `can_handle()` conservatively.
- Implement `compare()` and return `ConcordanceReport`.
- Keep `overall_concordance` in the 0-1 range.
- Add fixture files under `tests/fixtures/`.
- Add focused tests for perfect concordance and at least one discordant case.
- Document caveats, assumptions, and scoring behavior.

## Commit Style

Prefer small commits that complete one feature or fix. Run the full test suite
before committing.

