"""Strict base for authored protein-FASTA JSON documents."""

from pydantic import BaseModel, ConfigDict


class DocumentBase(BaseModel):
    """Forbid unknown document keys and post-validation mutation."""

    model_config = ConfigDict(extra="forbid", frozen=True)
