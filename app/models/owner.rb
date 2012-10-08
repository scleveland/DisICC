class Owner
  include DataMapper::Resource
  
  property :id,  Serial, :required => false, :key => true
  property :name, String, :required => false
  # property :first_name, String,  :required => false
  # property :last_name, String, :required => false
  # property :title, String, :required => false
  # property :email, String, :required => false
  # property :phone, String, :required => false

  has n, :items
  # has 1, :category
  # has 1, :owner
end