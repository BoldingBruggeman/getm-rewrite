class FortranObject:
    def get_arrays(self, names=(), default_dtype=float, dtypes={}, halo=None, **kwargs):
        for name in names:
            dtype = dtypes.get(name, default_dtype)
            data_ = self.get_array(name.encode('ascii'), dtype={float: 0, int: 1}[dtype], **kwargs)

            # Store field with and without halo as class attribute
            setattr(self, name + '_', data_)
            setattr(self, name, data_[tuple([slice(halo, -halo)] * data_.ndim)])