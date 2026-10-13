"""Explicit scope for the exact spherical V5.5 production branch."""
from dataclasses import dataclass


@dataclass(frozen=True)
class SphericalMatterScope:
    include_local_scalars: bool = True
    include_cosmos_scalar: bool = True
    include_dark_matter: bool = True
    include_baryons: bool = True
    include_radiation: bool = True
    include_shear: bool = False

    def validate(self) -> None:
        if self.include_shear:
            raise ValueError(
                "Spherical V5.5 cannot represent the archive's anisotropic "
                "Bianchi-I shear as an isotropic matter source."
            )
