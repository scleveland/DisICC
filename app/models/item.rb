class Item
  include DataMapper::Resource
  
  property :id,  Serial, :required => false, :key => true
  property :name, String,  :required => false
  property :description, String, :required => false
  property :url, String, :required => false
  property :contributor, String, :required => false
  property :review, Date, :required => false
  
  has n, :taggings
  has n, :topics, :through => :taggings
  belongs_to :category
    belongs_to :status
  belongs_to :owner
end