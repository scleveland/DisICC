class Topic
  include DataMapper::Resource
  
  property :id,  Serial, :required => false, :key => true
  property :name, String,  :required => false
  has n, :taggings
  has n,  :items, :through => :taggings
  validates_is_unique :name

end