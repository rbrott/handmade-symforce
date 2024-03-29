/** THIS IS AN AUTOMATICALLY GENERATED FILE.  DO NOT MODIFY
 * BY HAND!!
 *
 * Generated by lcm-gen
 *
 * From Source File: /Users/ryan/dev/symforce-zig-build/lcmtypes/symforce.lcm
 **/

#include <lcm/lcm_coretypes.h>

#ifndef __sym_linearization_offsets_t_hpp__
#define __sym_linearization_offsets_t_hpp__

#include <ostream>

namespace sym
{

class linearization_offsets_t
{
    public:
        // Offset of this key within the factor's state vector
        int32_t factor_offset;

        // Offset of this key within the whole problem's state vector
        int32_t combined_offset;

        // Tangent dimension of the key
        int32_t tangent_dim;

    public:
        linearization_offsets_t() = default;

        /**
         * Member constructor
         */
        inline linearization_offsets_t(
            const int32_t& factor_offset_arg,
            const int32_t& combined_offset_arg,
            const int32_t& tangent_dim_arg
        );

        /**
         * Encode a message into binary form.
         *
         * @param buf The output buffer.
         * @param offset Encoding starts at thie byte offset into @p buf.
         * @param maxlen Maximum number of bytes to write.  This should generally be
         *  equal to getEncodedSize().
         * @return The number of bytes encoded, or <0 on error.
         */
        inline __lcm_buffer_size encode(void *buf, __lcm_buffer_size offset, __lcm_buffer_size maxlen) const;

        /**
         * Check how many bytes are required to encode this message.
         */
        inline __lcm_buffer_size getEncodedSize() const;

        /**
         * Decode a message from binary form into this instance.
         *
         * @param buf The buffer containing the encoded message.
         * @param offset The byte offset into @p buf where the encoded message starts.
         * @param maxlen The maximum number of bytes to read while decoding.
         * @return The number of bytes decoded, or <0 if an error occured.
         */
        inline __lcm_buffer_size decode(const void *buf, __lcm_buffer_size offset, __lcm_buffer_size maxlen);

        /**
         * Retrieve the 64-bit fingerprint identifying the structure of the message.
         * Note that the fingerprint is the same for all instances of the same
         * message type, and is a fingerprint on the message type definition, not on
         * the message contents.
         */
        constexpr static int64_t getHash();

        using type_name_array_t = const char[24];

        inline static constexpr type_name_array_t* getTypeNameArrayPtr();

        /**
         * Returns "linearization_offsets_t"
         */
        inline static constexpr const char* getTypeName();

        /**
         * Returns "sym"
         */
        inline static constexpr const char * getPackageName();

        // LCM support functions. Users should not call these
        inline __lcm_buffer_size _encodeNoHash(void *buf, __lcm_buffer_size offset, __lcm_buffer_size maxlen) const;
        inline __lcm_buffer_size _getEncodedSizeNoHash() const;
        inline __lcm_buffer_size _decodeNoHash(const void *buf, __lcm_buffer_size offset, __lcm_buffer_size maxlen);
        constexpr static uint64_t _computeHash(const __lcm_hash_ptr *)
        {
            uint64_t hash = 0x320243c70940b2cfLL;
            return (hash<<1) + ((hash>>63)&1);
        }

        // Comparison operators.
        inline bool operator==(const linearization_offsets_t& other) const;
        inline bool operator!=(const linearization_offsets_t& other) const;

        // Ability to print to standard streams as well as the fmt library.
        friend std::ostream& operator<<(std::ostream& stream, const linearization_offsets_t& obj) {
#if defined(SKYMARSHAL_PRINTING_ENABLED)
            stream << "linearization_offsets_t(";
            stream << "factor_offset=" << obj.factor_offset << ", ";
            stream << "combined_offset=" << obj.combined_offset << ", ";
            stream << "tangent_dim=" << obj.tangent_dim;
            stream << ")";
#else
            stream << "<FORMATTING DISABLED>";
#endif
            return stream;
        }
};

linearization_offsets_t::linearization_offsets_t(
    const int32_t& factor_offset_arg,
    const int32_t& combined_offset_arg,
    const int32_t& tangent_dim_arg
) : factor_offset(factor_offset_arg),
    combined_offset(combined_offset_arg),
    tangent_dim(tangent_dim_arg) {}

__lcm_buffer_size linearization_offsets_t::encode(void *buf, __lcm_buffer_size offset, __lcm_buffer_size maxlen) const
{
    __lcm_buffer_size pos = 0, tlen;
    int64_t hash = (int64_t)getHash();

    tlen = __int64_t_encode_array(buf, offset + pos, maxlen - pos, &hash, 1);
    if(tlen < 0) return tlen; else pos += tlen;

    tlen = this->_encodeNoHash(buf, offset + pos, maxlen - pos);
    if (tlen < 0) return tlen; else pos += tlen;

    return pos;
}

__lcm_buffer_size linearization_offsets_t::decode(const void *buf, __lcm_buffer_size offset, __lcm_buffer_size maxlen)
{
    __lcm_buffer_size pos = 0, thislen;

    int64_t msg_hash;
    thislen = __int64_t_decode_array(buf, offset + pos, maxlen - pos, &msg_hash, 1);
    if (thislen < 0) return thislen; else pos += thislen;
    if (msg_hash != getHash()) return -1;

    thislen = this->_decodeNoHash(buf, offset + pos, maxlen - pos);
    if (thislen < 0) return thislen; else pos += thislen;

    return pos;
}

__lcm_buffer_size linearization_offsets_t::getEncodedSize() const
{
    return 8 + _getEncodedSizeNoHash();
}

constexpr int64_t linearization_offsets_t::getHash()
{
    return static_cast<int64_t>(_computeHash(NULL));
}

constexpr linearization_offsets_t::type_name_array_t* linearization_offsets_t::getTypeNameArrayPtr() {
    return &"linearization_offsets_t";
}

constexpr const char* linearization_offsets_t::getTypeName()
{
    return *linearization_offsets_t::getTypeNameArrayPtr();
}

constexpr const char * linearization_offsets_t::getPackageName()
{
    return "sym";
}

__lcm_buffer_size linearization_offsets_t::_encodeNoHash(void *buf, __lcm_buffer_size offset, __lcm_buffer_size maxlen) const
{
    __lcm_buffer_size pos = 0, tlen;

    tlen = __int32_t_encode_array(buf, offset + pos, maxlen - pos, &this->factor_offset, 1);
    if(tlen < 0) return tlen; else pos += tlen;

    tlen = __int32_t_encode_array(buf, offset + pos, maxlen - pos, &this->combined_offset, 1);
    if(tlen < 0) return tlen; else pos += tlen;

    tlen = __int32_t_encode_array(buf, offset + pos, maxlen - pos, &this->tangent_dim, 1);
    if(tlen < 0) return tlen; else pos += tlen;

    return pos;
}

__lcm_buffer_size linearization_offsets_t::_decodeNoHash(const void *buf, __lcm_buffer_size offset, __lcm_buffer_size maxlen)
{
    __lcm_buffer_size pos = 0, tlen;

    tlen = __int32_t_decode_array(buf, offset + pos, maxlen - pos, &this->factor_offset, 1);
    if(tlen < 0) return tlen; else pos += tlen;

    tlen = __int32_t_decode_array(buf, offset + pos, maxlen - pos, &this->combined_offset, 1);
    if(tlen < 0) return tlen; else pos += tlen;

    tlen = __int32_t_decode_array(buf, offset + pos, maxlen - pos, &this->tangent_dim, 1);
    if(tlen < 0) return tlen; else pos += tlen;

    return pos;
}

__lcm_buffer_size linearization_offsets_t::_getEncodedSizeNoHash() const
{
    __lcm_buffer_size enc_size = 0;
    enc_size += __int32_t_encoded_array_size(NULL, 1);
    enc_size += __int32_t_encoded_array_size(NULL, 1);
    enc_size += __int32_t_encoded_array_size(NULL, 1);
    return enc_size;
}

bool linearization_offsets_t::operator==(const linearization_offsets_t& other) const {
  return (
          (factor_offset==other.factor_offset) && 
          (combined_offset==other.combined_offset) && 
          (tangent_dim==other.tangent_dim));
}

bool linearization_offsets_t::operator!=(const linearization_offsets_t& other) const {
  return !(*this==other);
}

}

#endif