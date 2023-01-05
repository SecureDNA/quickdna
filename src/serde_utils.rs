/// Make "string-like" implementations of `serde::Deserialize` and
/// `serde::Serialize` for the given type, using that type's `Display` and
/// `FromStr` impls.
///
/// (This is used for example to serialize codons as `"ACG"` rather than
/// `["A","C","G"]`.)
macro_rules! impl_stringlike {
    ($type:ty) => {
        #[cfg(feature = "serde")]
        impl<'de> serde::Deserialize<'de> for $type {
            fn deserialize<D>(deserializer: D) -> Result<Self, D::Error>
            where
                D: serde::Deserializer<'de>,
            {
                struct Helper;
                impl<'de> serde::de::Visitor<'de> for Helper {
                    type Value = $type;

                    fn expecting(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
                        write!(formatter, "a string")
                    }

                    fn visit_str<E>(self, value: &str) -> Result<Self::Value, E>
                    where
                        E: serde::de::Error,
                    {
                        value
                            .parse::<Self::Value>()
                            .map_err(serde::de::Error::custom)
                    }
                }

                deserializer.deserialize_str(Helper)
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
