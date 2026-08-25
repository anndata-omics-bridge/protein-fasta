"""Shared strict base for authored protein-FASTA policy documents."""

from pydantic import BaseModel, ConfigDict


class PolicyDocument(BaseModel):
    """Forbid unknown policy keys and prevent post-validation mutation."""

    model_config = ConfigDict(extra="forbid", frozen=True)
