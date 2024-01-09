/** THIS IS AN AUTOMATICALLY GENERATED FILE.  DO NOT MODIFY
 * BY HAND!!
 *
 * Generated by lcm-gen
 *
 * From Source File: /Users/ryan/dev/symforce-zig-build/lcmtypes/symforce.lcm
 **/

#include <lcm/lcm_coretypes.h>

#ifndef __sym_index_entry_t_hpp__
#define __sym_index_entry_t_hpp__

#include <ostream>
#include "lcmtypes/sym/key_t.hpp"
#include "lcmtypes/sym/type_t.hpp"

namespace sym
{

class index_entry_t
{
    public:
        ::sym::key_t key;

        ::sym::type_t type;

        // Location within the storage or tangent vector, depending on context
        int32_t offset;

        // Size parameters
        int32_t storage_dim;

        int32_t tangent_dim;

    public:
        index_entry_t() = default;

        /**
         * Member constructor
         */
        inline index_entry_t(
            const ::sym::key_t& key_arg,
            const ::sym::type_t& type_arg,
            const int32_t& offset_arg,
            const int32_t& storage_dim_arg,
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
        inline static int64_t getHash();

        using type_name_array_t = const char[14];

        inline static constexpr type_name_array_t* getTypeNameArrayPtr();

        /**
         * Returns "index_entry_t"
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
#if !defined(SKYDIO_DISABLE_LCM_NO_INLINE)
        __attribute__((noinline))
#endif
        static uint64_t _computeHash(const __lcm_hash_ptr *p)
        {
            const __lcm_hash_ptr *fp;
            for(fp = p; fp != NULL; fp = fp->parent)
                if(fp->v == index_entry_t::getHash)
                    return 0;
            const __lcm_hash_ptr cp = { p, index_entry_t::getHash };

            uint64_t hash = 0xf55e0465bf8566a5LL +
                ::sym::key_t::_computeHash(&cp) +
         ::sym::type_t::_computeHash(&cp);

            return (hash<<1) + ((hash>>63)&1);
        }

        // Comparison operators.
        inline bool operator==(const index_entry_t& other) const;
        inline bool operator!=(const index_entry_t& other) const;

        // Ability to print to standard streams as well as the fmt library.
        friend std::ostream& operator<<(std::ostream& stream, const index_entry_t& obj) {
#if defined(SKYMARSHAL_PRINTING_ENABLED)
            stream << "index_entry_t(";
            stream << "key=" << obj.key << ", ";
            stream << "type=" << obj.type << ", ";
            stream << "offset=" << obj.offset << ", ";
            stream << "storage_dim=" << obj.storage_dim << ", ";
            stream << "tangent_dim=" << obj.tangent_dim;
            stream << ")";
#else
            stream << "<FORMATTING DISABLED>";
#endif
            return stream;
        }
};

index_entry_t::index_entry_t(
    const ::sym::key_t& key_arg,
    const ::sym::type_t& type_arg,
    const int32_t& offset_arg,
    const int32_t& storage_dim_arg,
    const int32_t& tangent_dim_arg
) : key(key_arg),
    type(type_arg),
    offset(offset_arg),
    storage_dim(storage_dim_arg),
    tangent_dim(tangent_dim_arg) {}

__lcm_buffer_size index_entry_t::encode(void *buf, __lcm_buffer_size offset, __lcm_buffer_size maxlen) const
{
    __lcm_buffer_size pos = 0, tlen;
    int64_t hash = (int64_t)getHash();

    tlen = __int64_t_encode_array(buf, offset + pos, maxlen - pos, &hash, 1);
    if(tlen < 0) return tlen; else pos += tlen;

    tlen = this->_encodeNoHash(buf, offset + pos, maxlen - pos);
    if (tlen < 0) return tlen; else pos += tlen;

    return pos;
}

__lcm_buffer_size index_entry_t::decode(const void *buf, __lcm_buffer_size offset, __lcm_buffer_size maxlen)
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

__lcm_buffer_size index_entry_t::getEncodedSize() const
{
    return 8 + _getEncodedSizeNoHash();
}

int64_t index_entry_t::getHash()
{
    static int64_t hash = _computeHash(NULL);
    return hash;
}

constexpr index_entry_t::type_name_array_t* index_entry_t::getTypeNameArrayPtr() {
    return &"index_entry_t";
}

constexpr const char* index_entry_t::getTypeName()
{
    return *index_entry_t::getTypeNameArrayPtr();
}

constexpr const char * index_entry_t::getPackageName()
{
    return "sym";
}

__lcm_buffer_size index_entry_t::_encodeNoHash(void *buf, __lcm_buffer_size offset, __lcm_buffer_size maxlen) const
{
    __lcm_buffer_size pos = 0, tlen;

    tlen = this->key._encodeNoHash(buf, offset + pos, maxlen - pos);
    if(tlen < 0) return tlen; else pos += tlen;

    tlen = this->type._encodeNoHash(buf, offset + pos, maxlen - pos);
    if(tlen < 0) return tlen; else pos += tlen;

    tlen = __int32_t_encode_array(buf, offset + pos, maxlen - pos, &this->offset, 1);
    if(tlen < 0) return tlen; else pos += tlen;

    tlen = __int32_t_encode_array(buf, offset + pos, maxlen - pos, &this->storage_dim, 1);
    if(tlen < 0) return tlen; else pos += tlen;

    tlen = __int32_t_encode_array(buf, offset + pos, maxlen - pos, &this->tangent_dim, 1);
    if(tlen < 0) return tlen; else pos += tlen;

    return pos;
}

__lcm_buffer_size index_entry_t::_decodeNoHash(const void *buf, __lcm_buffer_size offset, __lcm_buffer_size maxlen)
{
    __lcm_buffer_size pos = 0, tlen;

    tlen = this->key._decodeNoHash(buf, offset + pos, maxlen - pos);
    if(tlen < 0) return tlen; else pos += tlen;

    tlen = this->type._decodeNoHash(buf, offset + pos, maxlen - pos);
    if(tlen < 0) return tlen; else pos += tlen;

    tlen = __int32_t_decode_array(buf, offset + pos, maxlen - pos, &this->offset, 1);
    if(tlen < 0) return tlen; else pos += tlen;

    tlen = __int32_t_decode_array(buf, offset + pos, maxlen - pos, &this->storage_dim, 1);
    if(tlen < 0) return tlen; else pos += tlen;

    tlen = __int32_t_decode_array(buf, offset + pos, maxlen - pos, &this->tangent_dim, 1);
    if(tlen < 0) return tlen; else pos += tlen;

    return pos;
}

__lcm_buffer_size index_entry_t::_getEncodedSizeNoHash() const
{
    __lcm_buffer_size enc_size = 0;
    enc_size += this->key._getEncodedSizeNoHash();
    enc_size += this->type._getEncodedSizeNoHash();
    enc_size += __int32_t_encoded_array_size(NULL, 1);
    enc_size += __int32_t_encoded_array_size(NULL, 1);
    enc_size += __int32_t_encoded_array_size(NULL, 1);
    return enc_size;
}

bool index_entry_t::operator==(const index_entry_t& other) const {
  return (
          (key==other.key) && 
          (type==other.type) && 
          (offset==other.offset) && 
          (storage_dim==other.storage_dim) && 
          (tangent_dim==other.tangent_dim));
}

bool index_entry_t::operator!=(const index_entry_t& other) const {
  return !(*this==other);
}

}

#endif