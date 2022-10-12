# Wiki setup

## Install MkDocs
```bash
pip install mkdocs
```

Plugin to open images on the same page by *zooming* it:
```
pip install mkdocs-glightbox
```

To compile live on save
```
mkdocs serve
```

---

## Versioning

Consistency is TBD.
```bash
pip install mike

# make sure the default alias is 'latest'
mike deploy --push --update-aliases 0.0.0 latest
mike set-default --push latest
```

## Guidelines

Even though mkdocs is pure markdown, there are many extensions and plugins included. For the ones used here, refer to [mkdocs.yml](/mkdocs.yml). Furthermore, reference MkDocs [Setup](https://squidfunk.github.io/mkdocs-material/setup/changing-the-colors/) and [Reference](https://squidfunk.github.io/mkdocs-material/reference/).

## API Template

Refer to [API_template.md](/API_template.md) for the API structure.
