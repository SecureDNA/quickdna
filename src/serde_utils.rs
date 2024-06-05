// Copyright 2021-2024 SecureDNA Stiftung (SecureDNA Foundation) <licensing@securedna.org>
// SPDX-License-Identifier: MIT OR Apache-2.0

use std::{
    fmt::{self, Display},
    marker::PhantomData,
    str::FromStr,
};

/// A visitor object that uses T's FromStr impl to parse visited strings.
pub struct FromStrVisitor<T>(pub PhantomData<T>);

impl<'de, T: FromStr> serde::de::Visitor<'de> for FromStrVisitor<T>
where
    <T as FromStr>::Err: Display,
{
    type Value = T;

    fn expecting(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(formatter, "a string")
    }

    fn visit_str<E: serde::de::Error>(self, value: &str) -> Result<Self::Value, E> {
        value.parse().map_err(serde::de::Error::custom)
    }
}

/// Make "string-like" implementations of `serde::Deserialize` and
/// `serde::Serialize` for the given type, using that type's `Display` and
/// `FromStr` impls.
///
/// (This is used for example to serialize codons as `"ACG"` rather than
/// `["A","C","G"]`. We could equivalently use `serde_with::{DeserializeFromStr,
/// SerializeDisplay}`, but it's not worth the extra dependency.)
macro_rules! impl_stringlike {
    ($type:ty) => {
        #[cfg(feature = "serde")]
        impl<'de> serde::Deserialize<'de> for $type {
            fn deserialize<D>(deserializer: D) -> Result<Self, D::Error>
            where
                D: serde::Deserializer<'de>,
            {
                use crate::serde_utils::FromStrVisitor;
                deserializer.deserialize_str(FromStrVisitor(core::marker::PhantomData))
            }
        }

        #[cfg(feature = "serde")]
        impl serde::Serialize for $type {
            fn serialize<S>(&self, serializer: S) -> Result<S::Ok, S::Error>
            where
                S: serde::Serializer,
            {
                serializer.collect_str(self)
            }
        }
    };
}

pub(crate) use impl_stringlike;
