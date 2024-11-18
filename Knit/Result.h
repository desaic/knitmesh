#pragma once

#include <string>
#include <variant>

/// Less include-heavy version of std::errc. To convert, do
/// std::error_code(static_cast<int>(ErrcValue), std::generic_category())
enum class Errc {
  // success = 0, // always an error
  unknown = -1,
  address_family_not_supported = EAFNOSUPPORT,
  address_in_use = EADDRINUSE,
  address_not_available = EADDRNOTAVAIL,
  already_connected = EISCONN,
  argument_list_too_long = E2BIG,
  argument_out_of_domain = EDOM,
  bad_address = EFAULT,
  bad_file_descriptor = EBADF,
  bad_message = EBADMSG,
  broken_pipe = EPIPE,
  connection_aborted = ECONNABORTED,
  connection_already_in_progress = EALREADY,
  connection_refused = ECONNREFUSED,
  connection_reset = ECONNRESET,
  cross_device_link = EXDEV,
  destination_address_required = EDESTADDRREQ,
  device_or_resource_busy = EBUSY,
  directory_not_empty = ENOTEMPTY,
  executable_format_error = ENOEXEC,
  file_exists = EEXIST,
  file_too_large = EFBIG,
  filename_too_long = ENAMETOOLONG,
  function_not_supported = ENOSYS,
  host_unreachable = EHOSTUNREACH,
  identifier_removed = EIDRM,
  illegal_byte_sequence = EILSEQ,
  inappropriate_io_control_operation = ENOTTY,
  interrupted = EINTR,
  invalid_argument = EINVAL,
  invalid_seek = ESPIPE,
  io_error = EIO,
  is_a_directory = EISDIR,
  message_size = EMSGSIZE,
  network_down = ENETDOWN,
  network_reset = ENETRESET,
  network_unreachable = ENETUNREACH,
  no_buffer_space = ENOBUFS,
  no_child_process = ECHILD,
  no_link = ENOLINK,
  no_lock_available = ENOLCK,
  no_message = ENOMSG,
  no_protocol_option = ENOPROTOOPT,
  no_space_on_device = ENOSPC,
  no_stream_resources = ENOSR,
  no_such_device_or_address = ENXIO,
  no_such_device = ENODEV,
  no_such_file_or_directory = ENOENT,
  no_such_process = ESRCH,
  not_a_directory = ENOTDIR,
  not_a_socket = ENOTSOCK,
  not_a_stream = ENOSTR,
  not_connected = ENOTCONN,
  not_enough_memory = ENOMEM,
  not_supported = ENOTSUP,
  operation_canceled = ECANCELED,
  operation_in_progress = EINPROGRESS,
  operation_not_permitted = EPERM,
  operation_not_supported = EOPNOTSUPP,
  operation_would_block = EWOULDBLOCK,
  owner_dead = EOWNERDEAD,
  permission_denied = EACCES,
  protocol_error = EPROTO,
  protocol_not_supported = EPROTONOSUPPORT,
  read_only_file_system = EROFS,
  resource_deadlock_would_occur = EDEADLK,
  resource_unavailable_try_again = EAGAIN,
  result_out_of_range = ERANGE,
  state_not_recoverable = ENOTRECOVERABLE,
  stream_timeout = ETIME,
  text_file_busy = ETXTBSY,
  timed_out = ETIMEDOUT,
  too_many_files_open_in_system = ENFILE,
  too_many_files_open = EMFILE,
  too_many_links = EMLINK,
  too_many_symbolic_link_levels = ELOOP,
  value_too_large = EOVERFLOW,
  wrong_protocol_type = EPROTOTYPE
};

// maps to correct Errc type or unknown if out of range
Errc ErrcFromInt(int errc_int);

const char* ErrcToString(Errc);

struct Empty {
  auto operator<=>(const Empty&) const = default;
};

template <typename... Ts>
struct Overload : Ts... {
  using Ts::operator()...;
};
template <class... Ts>
Overload(Ts...) -> Overload<Ts...>;

class Error {
  std::variant<Errc, std::string> err_storage_;

 public:
  Error(Errc errc) : err_storage_(errc){};

  static Error from_errno();
  explicit Error(std::string owned_string) : err_storage_(std::move(owned_string)) {}

  const char* message() const {
    return std::visit(
        Overload{[](Errc errc) { return ErrcToString(errc); },
                 [](const std::string& owned_string) { return owned_string.c_str(); }},
        err_storage_);
  }

  Errc code() {
    if (std::holds_alternative<Errc>(err_storage_)) return std::get<Errc>(err_storage_);
    return Errc::unknown;
  }
};

template <typename ValueType>
  requires(!std::is_same_v<ValueType, Error>)
class [[nodiscard]] Result {
 private:
  using StorageType = std::variant<ValueType, Error>;
  StorageType storage_;

  template <typename OtherValueType>
    requires(!std::is_same_v<OtherValueType, Error>)
  friend class Result;

 public:
  Result()
    requires(std::same_as<ValueType, Empty>)
      : storage_(Empty{}) {}

  Result(Result&&) = default;
  Result& operator=(Result&&) = default;

  // Conditionally enable the copy constructor
  Result(const Result& other)
    requires(std::is_trivially_copyable_v<ValueType> || std::is_copy_constructible_v<ValueType>)
      : storage_(other.storage_) {}

  // Conditionally enable the copy assignment operator
  Result& operator=(const Result& other)
    requires(std::is_trivially_copyable_v<ValueType> || std::is_copy_constructible_v<ValueType>)
  {
    if (this != &other) {
      storage_ = other.storage_;
    }
    return *this;
  }

  template <typename U>
  Result(Result<U>&& other)
    requires(std::convertible_to<U, ValueType>)
      : storage_(std::visit(Overload{
                                [](U& v) { return StorageType{std::move(v)}; },
                                [](Error& e) { return StorageType{e}; },
                            },
                            other.storage_)) {}

  Result(ValueType value) : storage_(std::move(value)) {}
  Result(Errc errc) : storage_(Error(errc)) {}
  Result(Error error) : storage_(error) {}

  ValueType& value() { return std::get<0>(storage_); }

  ValueType const& value() const { return std::get<0>(storage_); }

  ValueType release_value() { return std::move(value()); }

  bool is_error() const { return storage_.index() == 1; }

  Error error() const { return std::get<1>(storage_); }

  ValueType release_value_but_fixme_should_propogate_errors() {
    // todo: add verify
    return release_value();
  }
};

template <>
class [[nodiscard]] Result<void> : public Result<Empty> {
 public:
  using Result<Empty>::Result;
};

#define CONCAT(a, b) CONCAT_INNER(a, b)
#define CONCAT_INNER(a, b) a##b
#define UNIQUE_NAME(base) CONCAT(base, __COUNTER__)

#define RETURN_IF_ERROR_IMPL(temporary_name, expr)                                       \
  do {                                                                                   \
    auto&& temporary_name = (expr);                                                      \
    static_assert(!std::is_lvalue_reference_v<decltype(temporary_name.release_value())>, \
                  "no references returned from fallible expressions");                   \
    if (temporary_name.is_error()) [[unlikely]]                                          \
      return temporary_name.error();                                                     \
  } while (0);

#define RETURN_IF_ERROR(expr) RETURN_IF_ERROR_IMPL(UNIQUE_NAME(_temporary_name), expr)

// @todo(mike) the below macros are awesome and only work in clang not msvc :(

// #define TRY_IMPL(temporary_name, expr)                                                   \
//   ({                                                                                     \
//     auto&& temporary_name = (expr);                                                      \
//     static_assert(!std::is_lvalue_reference_v<decltype(temporary_name.release_value())>, \
//                   "no references returned from fallible expressions");                   \
//     if (temporary_name.is_error()) [[unlikely]]                                          \
//       return temporary_name.error();                                                     \
//     temporary_name.release_value();                                                      \
//   })

// #define TRY(expr) TRY_IMPL(UNIQUE_NAME(_temporary_name), expr)
