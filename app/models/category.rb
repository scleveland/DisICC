class Category
  include DataMapper::Resource
  
  property :id,  Serial, :required => false, :key => true
  property :name, String,  :required => false

  has n,  :items



end