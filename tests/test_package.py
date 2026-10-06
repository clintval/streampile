import ast
import inspect
from importlib.resources import files

import streampile
from streampile import _native


def test_package_is_typed() -> None:
    assert files(streampile).joinpath("py.typed").is_file()


def test_the_native_stub_documents_each_object_as_its_docstring_does() -> None:
    stub = ast.parse(files(streampile).joinpath("_native.pyi").read_text())
    differences: list[str] = []

    def compare(node: ast.Module | ast.ClassDef, owner: object, prefix: str) -> None:
        for child in node.body:
            if not isinstance(child, ast.ClassDef | ast.FunctionDef):
                continue
            runtime = getattr(owner, child.name)
            documented = ast.get_docstring(child)
            if documented is not None and documented != inspect.cleandoc(runtime.__doc__ or ""):
                differences.append(prefix + child.name)
            if isinstance(child, ast.ClassDef):
                compare(child, runtime, f"{prefix}{child.name}.")

    compare(stub, _native, "")
    assert differences == []
