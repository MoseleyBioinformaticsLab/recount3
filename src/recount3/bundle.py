"""Resource bundles and concatenation helpers.

This module contains :class:`R3ResourceBundle`, a simple container for
workflows that operate on multiple :class:`R3Resource` objects.
"""

from __future__ import annotations

import dataclasses
from dataclasses import dataclass
from typing import Iterable, Tuple

import pandas as pd

from .bigwig import BigWigFile
from .errors import CompatibilityError
from .resource import R3Resource
from .types import CompatibilityMode


def _count_compat_keys(res: R3Resource) -> tuple[str, str]:
    """Return (family, feature_key) for a count resource."""
    rtype = getattr(res.description, "resource_type", None)
    if rtype == "count_files_gene_or_exon":
        gu = getattr(res.description, "genomic_unit", None) or ""
        family = "gene_or_exon"
        feature_key = f"{family}:{gu}"
        return family, feature_key

    if rtype == "count_files_junctions":
        jt = getattr(res.description, "junction_type", None) or ""
        jext = getattr(res.description, "junction_file_extension", None) or ""
        family = "junctions"
        feature_key = f"{family}:{jt}:{jext}"
        return family, feature_key

    raise ValueError(
        "Resource is not a recognized count-file type for stacking: "
        f"{rtype!r}"
    )


@dataclass(slots=True)
class R3ResourceBundle:
    """Container for a set of :class:`R3Resource` objects."""

    resources: list[R3Resource] = dataclasses.field(default_factory=list)

    def add(self, resource: R3Resource) -> None:
        """Add a resource to the bundle."""
        self.resources.append(resource)

    def extend(self, resources: Iterable[R3Resource]) -> None:
        """Extend the bundle with multiple resources."""
        self.resources.extend(resources)

    def load(self, *, strict: bool = True, force: bool = False) -> "R3ResourceBundle":
        """Load all resources and cache on each resource instance."""
        for res in self.resources:
            try:
                res.load(force=force)
            except Exception:
                if strict:
                    raise
                continue
        return self

    def iter_loaded(
        self,
        *,
        resource_type: str | None = None,
        autoload: bool = False,
    ) -> Iterable[tuple[R3Resource, object]]:
        """Yield (resource, object) pairs for resources with loaded data."""
        for res in self.resources:
            if resource_type and getattr(res.description, "resource_type", None) != resource_type:
                continue
            if not res.is_loaded():
                if autoload:
                    try:
                        res.load()
                    except Exception:
                        continue
                else:
                    continue
            obj = res.get_loaded()
            if obj is not None:
                yield res, obj

    def iter_bigwig(self, *, autoload: bool = True) -> Iterable[tuple[R3Resource, BigWigFile]]:
        """Yield (resource, BigWigFile) for BigWig resources."""
        for res, obj in self.iter_loaded(resource_type="bigwig_files", autoload=autoload):
            if isinstance(obj, BigWigFile):
                yield res, obj

    def get_loaded(self, *, resource_type: str | None = None, autoload: bool = False) -> list[object]:
        """Return loaded objects for resources in the bundle."""
        return [obj for _, obj in self.iter_loaded(resource_type=resource_type, autoload=autoload)]

    def filter(
        self,
        *,
        resource_type=None,
        organism=None,
        data_source=None,
        genomic_unit=None,
        project=None,
        sample=None,
        table_name=None,
        junction_type=None,
        predicate=None,
        invert: bool = False,
    ) -> "R3ResourceBundle":
        """Return a new bundle containing only resources that match criteria."""
        field_specs: dict[str, object] = {
            "resource_type": resource_type,
            "organism": organism,
            "data_source": data_source,
            "genomic_unit": genomic_unit,
            "project": project,
            "sample": sample,
            "table_name": table_name,
            "junction_type": junction_type,
        }
        field_specs = {k: v for k, v in field_specs.items() if v is not None}

        from .search import _match_spec  # local import to avoid cycles

        selected: list[R3Resource] = []
        for res in self.resources:
            desc = res.description
            fields_ok = all(_match_spec(getattr(desc, name, None), spec) for name, spec in field_specs.items())

            pred_ok = True
            if predicate is not None:
                try:
                    pred_ok = bool(predicate(res))
                except Exception:
                    pred_ok = False

            match = fields_ok and pred_ok
            if invert:
                match = not match

            if match:
                selected.append(res)

        return R3ResourceBundle(resources=selected)

    def only_counts(self) -> "R3ResourceBundle":
        """Return only gene/exon or junction count-file resources."""
        return self.filter(resource_type=("count_files_gene_or_exon", "count_files_junctions"))

    def only_metadata(self) -> "R3ResourceBundle":
        """Return only metadata-file resources."""
        return self.filter(resource_type="metadata_files")

    def exclude_metadata(self) -> "R3ResourceBundle":
        """Return a bundle with metadata-file resources removed."""
        return self.filter(resource_type="metadata_files", invert=True)

    def where(self, predicate) -> "R3ResourceBundle":
        """Predicate-based alias for ``filter(predicate=...)``."""
        return self.filter(predicate=predicate)

    # ------------------------------------------------------------------

    def stack_count_matrices(
        self,
        *,
        join: str = "inner",
        axis: int = 1,
        verify_integrity: bool = False,
        autoload: bool = True,
        compat: CompatibilityMode = "family",
    ):
        """Concatenate count matrices (gene/exon or junction) as DataFrames.

        Raises:
          CompatibilityError: If incompatible types are mixed.
          ImportError: If pandas is unavailable.
          TypeError: If a loaded object is not a DataFrame.
          ValueError: If no applicable resources are present.
        """
        if pd is None:
            raise ImportError("pandas is required for stack_count_matrices().")

        wanted = {"count_files_gene_or_exon", "count_files_junctions"}
        count_res = [r for r in self.resources if getattr(r.description, "resource_type", None) in wanted]
        if not count_res:
            raise ValueError("No count-file resources available to stack.")

        families: set[str] = set()
        features: set[str] = set()
        family_counts: dict[str, int] = {}

        for r in count_res:
            try:
                fam, feat = _count_compat_keys(r)
            except ValueError:
                continue
            families.add(fam)
            features.add(feat)
            family_counts[fam] = family_counts.get(fam, 0) + 1

        if compat == "family":
            if len(families) > 1:
                details = ", ".join(f"{k}={v}" for k, v in sorted(family_counts.items()))
                raise CompatibilityError(
                    "Incompatible count families selected for stacking. "
                    f"Found families: {sorted(families)} ({details}). "
                    "Stack gene/exon with gene/exon, and junctions with junctions. "
                    "Hint: filter first, e.g., "
                    'bundle.filter(resource_type="count_files_gene_or_exon") '
                    'or bundle.filter(resource_type="count_files_junctions").'
                )
        elif compat == "feature":
            if len(features) > 1:
                examples = ", ".join(sorted(features))
                raise CompatibilityError(
                    "Feature-level incompatibility detected. All inputs must share "
                    "the same feature key (e.g., gene vs exon; junction subtype). "
                    f"Distinct feature keys observed: {examples}. "
                    "Hint: filter by `genomic_unit` for gene/exon or by "
                    "`junction_type`/`junction_file_extension` for junctions."
                )
        else:
            raise ValueError(f"Unknown compat mode: {compat!r}")

        dfs: list[pd.DataFrame] = []
        for res, obj in self.iter_loaded(autoload=autoload):
            rtype = getattr(res.description, "resource_type", None)
            if rtype not in wanted:
                continue
            if not isinstance(obj, pd.DataFrame):
                raise TypeError(f"Loaded object for {res.url} is not a DataFrame.")
            dfs.append(obj)

        if not dfs:
            raise ValueError(
                "No loaded count matrices found. "
                "Try autoload=True or call bundle.load() first."
            )

        return pd.concat(dfs, axis=axis, join=join, verify_integrity=verify_integrity)  # type: ignore
