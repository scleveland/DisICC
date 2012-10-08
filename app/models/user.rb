class User

  include DataMapper::Resource

  property :id, Serial

  property :first_nam, String
  property :last_name, String
  property :password, String

end
