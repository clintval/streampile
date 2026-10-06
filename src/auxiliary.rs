use crate::error::{Result, invalid_data};

/// The value of a record's auxiliary field, borrowed from the record without copying.
#[derive(Clone, Copy, Debug, PartialEq)]
pub enum AuxValue<'a> {
    /// A printable character, type `A`.
    Character(u8),
    /// An integer of any width, types `c`, `C`, `s`, `S`, `i`, and `I`.
    Integer(i64),
    /// A single-precision float, type `f`.
    Float(f32),
    /// A string, type `Z`, without its terminating NUL.
    String(&'a [u8]),
    /// A hex string, type `H`, without its terminating NUL.
    Hex(&'a [u8]),
    /// A numeric array, type `B`.
    Array(AuxArray<'a>),
}

impl AuxValue<'_> {
    /// The element at an index of a string or an array, such as a per-base tag at a query
    /// offset, or `None` for a scalar or an index past the end.
    pub fn get(&self, index: usize) -> Option<AuxElement> {
        match self {
            AuxValue::String(bytes) | AuxValue::Hex(bytes) => {
                bytes.get(index).copied().map(AuxElement::Byte)
            }
            AuxValue::Array(array) => array.get(index),
            _ => None,
        }
    }
}

/// One element of a string or array auxiliary field.
#[derive(Clone, Copy, Debug, PartialEq)]
pub enum AuxElement {
    /// A byte of a string.
    Byte(u8),
    /// An element of an integer array.
    Integer(i64),
    /// An element of a float array.
    Float(f32),
}

/// The element type of a numeric array auxiliary field.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum ArraySubtype {
    /// `c`
    Int8,
    /// `C`
    UInt8,
    /// `s`
    Int16,
    /// `S`
    UInt16,
    /// `i`
    Int32,
    /// `I`
    UInt32,
    /// `f`
    Float,
}

impl ArraySubtype {
    fn from_byte(byte: u8) -> Result<Self> {
        Ok(match byte {
            b'c' => ArraySubtype::Int8,
            b'C' => ArraySubtype::UInt8,
            b's' => ArraySubtype::Int16,
            b'S' => ArraySubtype::UInt16,
            b'i' => ArraySubtype::Int32,
            b'I' => ArraySubtype::UInt32,
            b'f' => ArraySubtype::Float,
            other => {
                return Err(invalid_data(format!(
                    "invalid array subtype {}",
                    other.escape_ascii()
                ))
                .into());
            }
        })
    }

    /// The width of one element in bytes.
    pub fn width(self) -> usize {
        match self {
            ArraySubtype::Int8 | ArraySubtype::UInt8 => 1,
            ArraySubtype::Int16 | ArraySubtype::UInt16 => 2,
            ArraySubtype::Int32 | ArraySubtype::UInt32 | ArraySubtype::Float => 4,
        }
    }
}

/// A numeric array auxiliary field, read element by element from the record's bytes.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct AuxArray<'a> {
    subtype: ArraySubtype,
    bytes: &'a [u8],
}

impl<'a> AuxArray<'a> {
    /// The element type.
    pub fn subtype(&self) -> ArraySubtype {
        self.subtype
    }

    /// The elements as little-endian bytes, as stored.
    pub fn as_bytes(&self) -> &'a [u8] {
        self.bytes
    }

    /// The number of elements.
    pub fn len(&self) -> usize {
        self.bytes.len() / self.subtype.width()
    }

    /// Whether the array has no elements.
    pub fn is_empty(&self) -> bool {
        self.bytes.is_empty()
    }

    /// The element at an index, or `None` past the end.
    pub fn get(&self, index: usize) -> Option<AuxElement> {
        let width = self.subtype.width();
        let start = index.checked_mul(width)?;
        let bytes = self.bytes.get(start..start.checked_add(width)?)?;
        Some(element(self.subtype, bytes))
    }

    /// Every element, in order.
    pub fn iter(&self) -> impl ExactSizeIterator<Item = AuxElement> + 'a {
        let subtype = self.subtype;
        self.bytes
            .chunks_exact(subtype.width())
            .map(move |bytes| element(subtype, bytes))
    }
}

fn element(subtype: ArraySubtype, bytes: &[u8]) -> AuxElement {
    match subtype {
        ArraySubtype::Int8 => AuxElement::Integer(i64::from(bytes[0] as i8)),
        ArraySubtype::UInt8 => AuxElement::Integer(i64::from(bytes[0])),
        ArraySubtype::Int16 => {
            AuxElement::Integer(i64::from(i16::from_le_bytes([bytes[0], bytes[1]])))
        }
        ArraySubtype::UInt16 => {
            AuxElement::Integer(i64::from(u16::from_le_bytes([bytes[0], bytes[1]])))
        }
        ArraySubtype::Int32 => AuxElement::Integer(i64::from(i32::from_le_bytes(word(bytes)))),
        ArraySubtype::UInt32 => AuxElement::Integer(i64::from(u32::from_le_bytes(word(bytes)))),
        ArraySubtype::Float => AuxElement::Float(f32::from_le_bytes(word(bytes))),
    }
}

fn word(bytes: &[u8]) -> [u8; 4] {
    [bytes[0], bytes[1], bytes[2], bytes[3]]
}

/// The byte range of one field's type and value within a record's auxiliary data.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub(crate) struct Field {
    pub start: u32,
    pub end: u32,
}

impl Field {
    /// The value of the field, decoded again from data it was found in.
    pub fn value(self, data: &[u8]) -> Result<AuxValue<'_>> {
        let bytes = data
            .get(self.start as usize..self.end as usize)
            .ok_or_else(|| invalid_data("an auxiliary field is out of range"))?;
        Ok(decode(bytes)?.0)
    }
}

/// Walks every field of a record's auxiliary data, calling `visit` with each tag and the range
/// of its type and value.
pub(crate) fn walk(data: &[u8], mut visit: impl FnMut([u8; 2], Field)) -> Result<()> {
    for field in fields(data) {
        let (tag, field) = field?;
        visit(tag, field);
    }
    Ok(())
}

/// Finds the value of the first field with a tag in a record's auxiliary data, reading no
/// further than it.
pub(crate) fn find(data: &[u8], tag: [u8; 2]) -> Result<Option<AuxValue<'_>>> {
    for field in fields(data) {
        let (found, field) = field?;
        if found == tag {
            return field.value(data).map(Some);
        }
    }
    Ok(None)
}

/// Each field of a record's auxiliary data with its tag, in order, ending at the first error.
fn fields(data: &[u8]) -> impl Iterator<Item = Result<([u8; 2], Field)>> + '_ {
    let mut at = 0;
    std::iter::from_fn(move || {
        if at >= data.len() {
            return None;
        }
        let field = next_field(data, at);
        at = field
            .as_ref()
            .map_or(data.len(), |(_, field)| field.end as usize);
        Some(field)
    })
}

fn next_field(data: &[u8], at: usize) -> Result<([u8; 2], Field)> {
    let tag = data
        .get(at..at + 2)
        .ok_or_else(|| invalid_data("an auxiliary field is truncated"))?;
    let start = at + 2;
    let (_, length) = decode(&data[start.min(data.len())..])?;
    let field = Field {
        start: u32::try_from(start).map_err(|_| invalid_data("auxiliary data is too long"))?,
        end: u32::try_from(start + length)
            .map_err(|_| invalid_data("auxiliary data is too long"))?,
    };
    Ok(([tag[0], tag[1]], field))
}

/// Decodes a type byte and the value after it, and returns the value and the bytes it took.
fn decode(src: &[u8]) -> Result<(AuxValue<'_>, usize)> {
    let truncated = || invalid_data("an auxiliary field is truncated");
    let (&kind, value) = src.split_first().ok_or_else(truncated)?;
    let fixed = |width: usize| value.get(..width).ok_or_else(truncated);
    Ok(match kind {
        b'A' => (AuxValue::Character(fixed(1)?[0]), 2),
        b'c' => (AuxValue::Integer(i64::from(fixed(1)?[0] as i8)), 2),
        b'C' => (AuxValue::Integer(i64::from(fixed(1)?[0])), 2),
        b's' => {
            let bytes = fixed(2)?;
            (
                AuxValue::Integer(i64::from(i16::from_le_bytes([bytes[0], bytes[1]]))),
                3,
            )
        }
        b'S' => {
            let bytes = fixed(2)?;
            (
                AuxValue::Integer(i64::from(u16::from_le_bytes([bytes[0], bytes[1]]))),
                3,
            )
        }
        b'i' => (
            AuxValue::Integer(i64::from(i32::from_le_bytes(word(fixed(4)?)))),
            5,
        ),
        b'I' => (
            AuxValue::Integer(i64::from(u32::from_le_bytes(word(fixed(4)?)))),
            5,
        ),
        b'f' => (AuxValue::Float(f32::from_le_bytes(word(fixed(4)?))), 5),
        b'Z' | b'H' => {
            let length = value
                .iter()
                .position(|&byte| byte == 0)
                .ok_or_else(truncated)?;
            let bytes = &value[..length];
            let decoded = if kind == b'Z' {
                AuxValue::String(bytes)
            } else {
                AuxValue::Hex(bytes)
            };
            (decoded, length + 2)
        }
        b'B' => {
            let (&subtype, rest) = value.split_first().ok_or_else(truncated)?;
            let subtype = ArraySubtype::from_byte(subtype)?;
            let count = rest.get(..4).ok_or_else(truncated)?;
            let count = u32::from_le_bytes(word(count)) as usize;
            let length = count.checked_mul(subtype.width()).ok_or_else(truncated)?;
            let end = length.checked_add(4).ok_or_else(truncated)?;
            let bytes = rest.get(4..end).ok_or_else(truncated)?;
            (AuxValue::Array(AuxArray { subtype, bytes }), end + 2)
        }
        other => {
            return Err(invalid_data(format!(
                "invalid auxiliary field type {}",
                other.escape_ascii()
            ))
            .into());
        }
    })
}
